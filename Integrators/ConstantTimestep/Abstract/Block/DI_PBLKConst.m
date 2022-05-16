
% ======================================================================================================================
%  Parallel Diagonally Implicit Block Method : y^{(n+1)} = A * y^{(n)} + r * B * y^{(n+1)} + C * y^{(n+1)} + r * D * f^{[n+1]}
%  with optional output point
%       y_out = a_out * y^{(n)} + r * b_out * y^{(n+1)} + c_out * y^{(n+1)} + r * d_out * f^{[n+1]} + r * e_out * f_out
%
%   A parallel diagonally implicit block method must satisfy the following conditions:
%       1. The matrix C must be the zero matrix.
%       2. The matrix D must be a diagonal matrix.
% ======================================================================================================================

classdef DI_PBLKConst < DI_BLKConst
    
    methods
        
        function this = DI_PBLKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@DI_BLKConst(options);
        end
        
    end
    
    methods
        
        function [t_out, y_out] = solve(this, problem) % == redefine solve to allow for spmd ===========================
            if(isempty(this.num_timesteps))
                error('Integrator: number of timesteps has not been set.');
            end
            % -- set stepsize ------------------------------------------------------------------------------------------
            this.setStepsize(problem);
            % -- obtain initial conditions -----------------------------------------------------------------------------
            [t_out, y_out, clean_exit] = this.initialConditions(problem);
            if(~clean_exit)
                t_out = NaN;
                y_out = NaN;
                warning('Integrator: Emergency exit during initial condition computation.');
                return;
            end
            
            num_proc  = this.getNumProcessors();
            
            % spmd "break-in" : compute reference solution (or this will be automatically done when problem is passed to
            % workers which artificially slows down time. Also open empty pool in case one is not already open.
            temp = problem.initial_condition;
            spmd(num_proc)
            end
            
            % -- solve problem -----------------------------------------------------------------------------------------
            start_step_time   = tic;
            num_steps = this.num_timesteps;
            
            spmd(num_proc)
                [struct_out, y_out] = this.initStepStruct(t_out, y_out, problem);
                for i = 1 : num_steps
                    % --> step
                    t_in      = t_out;
                    struct_in = struct_out;
                    y_in      = y_out;                
                    [t_out, y_out, struct_out] = this.step(t_in, y_in, struct_in, problem);
                    % --> emergency exit conditions
                    if(labindex == num_proc)
                        emergency_exit_flag = any(isinf(y_in(:))) || any(isnan(y_in(:)));
                        labSend(emergency_exit_flag, 1:num_proc-1, 2)
                        if(emergency_exit_flag)
                            t_out = NaN;
                            y_out = NaN;
                            warning('Integrator: Emergency exit at step %i', i);
                            break;
                        end
                    else
                        emergency_exit_flag = labReceive(num_proc, 2);
                        if(emergency_exit_flag)
                            t_out = NaN;
                            y_out = NaN;
                            break;
                        end                        
                    end
                end
            end
            [t_out, y_out] = this.userOutput(t_in{num_proc}, y_in{num_proc}, struct_in{num_proc}, t_out{num_proc}, y_out{num_proc}, struct_out{num_proc}, problem);
            
            this.step_stats.recordStep(toc(start_step_time));
            % -- output solution if only one output argument -----------------------------------------------------------
            if(nargout <= 1)
                t_out = y_out;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, ~, y_in, problem)
            
            num_proc = this.getNumProcessors();
            if(labindex == num_proc)
                if(~this.parallel_initial_guess) % ensure that parallel_initial_guess is true
                    this.parallel_initial_guess = true;
                    warning('parallel_initial_guess set to true. This flag must always be true for parallel method.')
                end
                [step_struct, y_in] = initStepStruct@DI_BLKConst(this, [], y_in, problem);
                step_struct.conj_inds  = find(this.conjugate_outputs(:).' ~= 0); %ensure row vector so that it can be used as for loop indices
            else
                step_struct = struct();
            end
            % -- additional parameters for all threads
            step_struct.num_proc = num_proc;
            step_struct.solve_inds = find(this.conjugate_outputs == 0);
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            num_proc   = step_struct.num_proc;
            solve_inds = step_struct.solve_inds;
            
            if(labindex == num_proc) % leader processor : compute nonlinear system vectors
                
                % -- read data from struct -----------------------------------------------------------------------------
                h        = this.h;
                q        = step_struct.q;
                AT       = step_struct.AT;
                BT       = step_struct.BT;
                AT_extrp = step_struct.AT_extrp;
                BT_extrp = step_struct.BT_extrp;
                noi      = step_struct.noi;
                req_RHS_flags = step_struct.req_RHS_F;
                
                t_out = t_in + h;
                y_out = zeros(size(y_in));
                b_j   = zeros(size(y_in, 1), 1);
                num_solves = num_proc - 1;
                
                % -- create empty time record for current step ---------------------------------------------------------
                for j = 1 : q
                    this.step_stats.setIndex(j);
                    this.step_stats.recordStep(0);
                end
                
                % -- compute and send information for nonlinear system y = b + c * F(y) --------------------------------
                for j = 1 : num_solves
                    
                    output_index = solve_inds(j);
                    this.step_stats.startTimer(output_index);
                    
                    if(problem.real_valued == false || this.conjugate_outputs(output_index) == 0)
                        this.nonlinear_solver.stats.setIndex(output_index);
                        this.nonlinear_solver.linear_solver.stats.setIndex(output_index);
                        
                        % -- construct vector b ------------------------------------------------------------------------
                        %  b_j = y_in * AT(:, j) + h * f_in * BT(:, j) + y_out * CT(:, j) + f_out * DT(:, 1:j-1)
                        % ----------------------------------------------------------------------------------------------
                        b_j(:) = 0;
                        if(~isempty(this.non_zero_A_indices{output_index}))
                            b_j = b_j + y_in(:, this.non_zero_A_indices{output_index}) * AT(this.non_zero_A_indices{output_index}, output_index);
                        end
                        if(~isempty(this.non_zero_B_indices{output_index}))
                            b_j = b_j + h * (step_struct.F_in(:, this.non_zero_B_indices{output_index}) * BT(this.non_zero_B_indices{output_index}, output_index));
                        end
                        % Note: C = 0 and tril(D) = 0 since method is parallel; therfore they can be skipped.
                        
                        % -- set constant c ----------------------------------------------------------------------------
                        c_j = h * this.D(output_index, output_index);
                        
                        % -- choose guess for nonlinear system ---------------------------------------------------------
                        if(this.extrapolate_initial_guess)
                            y_guess = zeros(size(y_in, 1), 1);
                            if(~isempty(this.A_extrapolate))
                                y_guess = y_guess + y_in * AT_extrp(:, output_index);
                            end
                            if(~isempty(this.B_extrapolate))
                                y_guess = y_guess + h * step_struct.F_in * BT_extrp(:, output_index);
                            end
                        else % choose nearest value as guess -----------------------------------------------------------
                            if(noi(output_index) <= q)
                                y_guess = y_in(:, noi(output_index));
                            else
                                y_guess = y_out(:, noi(output_index) - q);
                            end
                        end
                        
                        eval_f_flag = req_RHS_flags(output_index);
                        labSend({b_j, c_j, y_guess, eval_f_flag}, j, 0);
                        
                    end
                    
                end
                
                % -- recieve solutions to nonlinear systems y = b + c * F(y) -------------------------------------------
                clean_exit = true;
                for j = 1 : num_solves
                    [data, srcWkrIdx] = labReceive('any', 1);
                    output_index = solve_inds(srcWkrIdx);
                    y_out(:, output_index) = data{1};
                    if(req_RHS_flags(output_index))
                        step_struct.F_out(:, output_index) = data{2};
                    end
                    clean_exit = clean_exit && data{3};
                end
                
                if(clean_exit)
                    % -- set conjugate outputs -------------------------------------------------------------------------
                    for conj_ind = step_struct.conj_inds
                        y_out(:, conj_ind) = conj(y_out(:, this.conjugate_outputs(conj_ind)));
                        step_struct.F_out(:, conj_ind) = conj(step_struct.F_out(:, this.conjugate_outputs(conj_ind)));
                    end

                    % -- update F_in for next step ---------------------------------------------------------------------
                    step_struct.F_in = step_struct.F_out;
                else
                    y_out = NaN;
                    t_out = NaN;
                end
                
            else % worker processor : solve nonlinear system
                
                output_index = solve_inds(labindex);
                
                % -- recieve information for nonlinear system y = b + c * F(y) -------------------------------------
                data = labReceive(num_proc, 0);
                b  = data{1};
                c  = data{2};
                y0 = data{3};
                eval_f_flag = data{4};
                
                % -- solve nonlinear system ------------------------------------------------------------------------
                if(this.real_valued_outputs(output_index)) % real valued
                    if(c ~= 0)
                        [y, clean_exit] = this.nonlinear_solver.solveBC(problem, real(b), real(c), real(y0));
                    else
                        y = real(b);
                    end
                else % complex valued
                    if(c ~= 0)
                        [y, clean_exit] = this.nonlinear_solver.solveBC(problem, b, c, y0);
                    else
                        y = b;
                    end
                end
                
                % -- evaluate RHS ----------------------------------------------------------------------------------
                if(eval_f_flag)
                    this.rhs_stats.startTimer();
                    if(this.eval_RHS)
                        f = problem.RHS(y);
                        this.rhs_stats.recordRHSEval();
                    else
                        if(this.real_valued_outputs(output_index)) % real valued
                            f = real(y - b)/real(c);
                        else
                            f = (y - b)/c;
                        end
                    end
                    this.rhs_stats.recordRHSEval();
                else
                    f = [];
                end
                
                % -- send solution ---------------------------------------------------------------------------------
                labSend({y,f,clean_exit}, num_proc, 1);
                y_out = [];
                t_out = [];
                
            end
            
        end
    end
    
    methods (Access = protected)
        
        function num_proc = getNumProcessors(this)
            num_proc = 1 + sum(this.conjugate_outputs == 0);
        end
        
        function verifyCefficientMatrices(this)
            if(any(this.C))
                error('Coefficient Matrix C must be empty for parallel block method.');
            end
            if(~isdiag(this.D))
                error('Coefficient Matrix D must be diagonal for parallel block method.');
            end
        end
        
    end
    
end