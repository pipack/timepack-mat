
% ======================================================================================================================
%
%  Diagonally Implicit Block Method:
%       y^{(n+1)} = A * y^{(n)} + r * B * y^{(n+1)} + C * y^{(n+1)} + r * D * f^{[n+1]}
%  with optional output point 
%       y_out = a_out * y^{(n)} + r * b_out * f^{(n)} + c_out * y^{(n+1)} + r * d_out * f^{[n+1]} + r * e_out * f_ou 
%   
%  A diagonally implicit block method must satisfy the following conditions:
%   1. The matrix D must be lower triagular
%
% ======================================================================================================================

classdef DI_BLKConst < BLKConst & ImplicitIntegratorConst
    
    methods
        
        function this = DI_BLKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@BLKConst(options);
            this = this@ImplicitIntegratorConst(options);
            this.verifyCefficientMatrices();
            this.setNonZeroMatrixIndices();
        end
        
        function setProperties(this, prop_struct)
            setProperties@BLKConst(this, prop_struct);
            
            % -- recompute fields if coefficient matrices are redefined ------------------------------------------------
            fields = fieldnames(prop_struct);
            if(contains(fields, {'A', 'B', 'C', 'D'})) 
                this.verifyCefficientMatrices();
                this.setNonZeroMatrixIndices();
            end
        end
        
    end
    
    methods (Access = protected)
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            % -- read data from struct ---------------------------------------------------------------------------------          
            h     = this.h;
            q     = step_struct.q;
            AT    = step_struct.AT;
            BT    = step_struct.BT;
            CT    = step_struct.CT;
            DT    = step_struct.DT;
            AT_extrp = step_struct.AT_extrp;
            BT_extrp = step_struct.BT_extrp;
            noi   = step_struct.noi;
            req_RHS_flags = step_struct.req_RHS_F;
            
            t_out = t_in + h;
            y_out = zeros(size(y_in));
            b_j   = zeros(size(y_in, 1), 1);
            
            % -- compute outputs ---------------------------------------------------------------------------------------
            for j = 1 : q
            	
                this.step_stats.startTimer(j);
                
                if(problem.real_valued == false || this.conjugate_outputs(j) == 0)
                    this.nonlinear_solver.stats.setIndex(j);
                    this.nonlinear_solver.linear_solver.stats.setIndex(j);
                        
                    % -- b_j = y_in * AT(:, j) + h * f_in * BT(:, j) + y_out * CT(:, j) + f_out * DT(:, 1:j-1) ---------
                    b_j(:) = 0;    
                    if(~isempty(this.non_zero_A_indices{j}))
                        b_j = b_j + y_in(:, this.non_zero_A_indices{j}) * AT(this.non_zero_A_indices{j}, j);
                    end                    
                    if(~isempty(this.non_zero_B_indices{j}))
                        b_j = b_j + h * (step_struct.F_in(:, this.non_zero_B_indices{j}) * BT(this.non_zero_B_indices{j}, j));
                    end
                    if(~isempty(this.non_zero_C_indices{j}))
                        b_j = b_j + y_out(:, this.non_zero_C_indices{j}) * CT(this.non_zero_C_indices{j}, j);
                    end
                    if(~isempty(this.non_zero_D_indices{j}))
                        b_j = b_j + h * (step_struct.F_out(:, this.non_zero_D_indices{j}))  * DT(this.non_zero_D_indices{j}, j);
                    end
                    
                    c_j = h * this.D(j, j);
                    
                    if(c_j ~= 0) % -- solve nonlinear system y_j = b_j + c_j * f(y_j) ----------------------------------
                        
                        % -- choose guess for nonlinear system ---------------------------------------------------------
                        if(this.extrapolate_initial_guess)
                            y_guess = zeros(size(y_in, 1), 1);
                            if(~isempty(this.A_extrapolate))
                                y_guess = y_guess + y_in * AT_extrp(:, j);
                            end
                            if(~isempty(this.B_extrapolate))
                                y_guess = y_guess + h * step_struct.F_in * BT_extrp(:, j);
                            end
                        else % choose nearest value as guess
                            if(noi(j) <= q)
                                y_guess = y_in(:, noi(j));
                            else
                                y_guess = y_out(:, noi(j) - q);
                            end
                        end
                        
                        if(problem.real_valued && this.real_valued_outputs(j)) % -- output is real valued time point -------------
                        	[y_out(:,j), clean_exit] = this.nonlinear_solver.solveBC(problem, real(b_j), real(c_j), real(y_guess));
                        else % -- output is complex valued -------------------------------------------------------------
                        	[y_out(:,j), clean_exit] = this.nonlinear_solver.solveBC(problem, b_j, c_j, y_guess);
                        end

                        emergency_exit = ((~clean_exit) || any(isinf(y_out(:,j))) || any(isnan(y_out(:,j))));
                        if(emergency_exit)
                            y_out = NaN;
                            return;
                        end

                    else % -- system is explicit -----------------------------------------------------------------------
                        if(problem.real_valued && this.real_valued_outputs(j))
                            y_out(:,j) = real(b_j);
                        else
                            y_out(:,j) = b_j;
                        end
                    end
                
                else
                    conj_j = this.conjugate_outputs(j);
                    y_out(:, j) = conj(y_out(:, conj_j));
                end
                    
                % -- eval RHS ------------------------------------------------------------------------------------------
                this.rhs_stats.startTimer(j);
                if(req_RHS_flags(j))
                    conj_j = this.conjugate_outputs(j);
                    if(conj_j ~= 0 && problem.real_valued) % -- check for conjugate input ------------------------------
                        step_struct.F_out(:, j) = conj(step_struct.F_out(:, conj_j));
                    else
                        if(this.eval_RHS || c_j == 0)
                            step_struct.F_out(:, j) = problem.RHS(y_out(:, j)); % -- eval directly ---------------------
                        else % -- compute via re-arrangement -----------------------------------------------------------
                            if(problem.real_valued && this.real_valued_outputs(j))
                                step_struct.F_out(:, j) = real(y_out(:,j) - b_j) / real(c_j); 
                            else
                                step_struct.F_out(:, j) = (y_out(:,j) - b_j) / c_j; 
                            end
                        end
                    end
                    
                end
                this.rhs_stats.recordRHSEval();
                this.step_stats.recordStep();
            end
            
            if(final_step && ~isempty(this.a_out))
                [t_out, y_out] = this.finalStep(t_in, y_in, y_out, step_struct, problem);
            end
            
            % -- update F_in for next step -----------------------------------------------------------------------------
            step_struct.F_in = step_struct.F_out;
            
        end
        
        function [t_final, y_final] = finalStep(this, t_in, y_in, y_out, step_struct, problem)
            
            flag = ~isempty(this.a_out) && ~isempty(this.b_out) && ~isempty(this.c_out) && ~isempty(this.d_out) && ~isempty(this.e_out);
            if(flag)
                
                % -- read data from struct ---------------------------------------------------------------------------------
                h = this.h;
                q = step_struct.q;
                req_RHS_flags = step_struct.req_RHS_F;
                
                % -- check if F input requires additional Evaluations ------------------------------------------------------
                uncomputed_F_in  = setdiff(find(this.b_out), find(req_RHS_flags));
                for j = uncomputed_F_in(:)'
                    this.rhs_stats.startTimer(j);
                    step_struct.F_in(:, j) = problem.RHS(y_in(:, j));
                    this.rhs_stats.recordRHSEval();
                end
                % -- check if F output requires additional Evaluations -----------------------------------------------------
                uncomputed_F_out = setdiff(find(this.d_out), find(req_RHS_flags));
                for j = uncomputed_F_out(:)'
                    this.rhs_stats.startTimer(j);
                    step_struct.F_out(:, j) = problem.RHS(y_out(:, j));
                    this.rhs_stats.recordRHSEval();
                end
                
                % -- form nonlinear equation y = b + c * f(y) --------------------------------------------------------------
                b = y_in * this.a_out(:) + h * (step_struct.F_in * this.b_out(:)) + ...
                    y_out * this.c_out(:) + h * (step_struct.F_out * this.d_out(:));
                c = h * this.e_out;
                
                if(c ~= 0) % -- solve nonlinear system ---------------------------------------------------------------------
                    
                    % -- choose guess for nonlinear system -----------------------------------------------------------------
                    nearest_index = this.nearestPoint([this.starting_times(:); this.starting_times(:) + this.h], this.z_out + this.h);
                    if(nearest_index <= q)
                        y_guess = y_in(:, nearest_index);
                    else
                        y_guess = y_out(:, nearest_index - q);
                    end
                    % -- solve system --------------------------------------------------------------------------------------
                    [y_final, clean_exit] = this.nonlinear_solver.solveBC(problem, b, c, y_guess);
                    if(imag(this.z_out) == 0 && problem.real_valued)
                        y_final = real(y_final);
                    end
                    emergency_exit = ((~clean_exit) || any(isinf(y_final)) || any(isnan(y_final)));
                    if(emergency_exit)
                        t_final = NaN;
                        y_final = NaN;
                        return;
                    end
                else % -- system is explicit -------------------------------------------------------------------------------
                    if(imag(this.z_out) == 0 && problem.real_valued)
                        y_final = real(b);
                    else
                        y_final = b;
                    end
                end
                
                t_final = problem.tspan(end) + h * this.z_out; % output value
                
            else
                t_final = t_in + this.h;
                y_final = y_in;
            end
            
        end
        
    end
    
    methods (Access = protected)
        
        function verifyCefficientMatrices(this)
            isstril = @(A) istril(A) && ~any(diag(A)); % tests if A is strictly lower triangular
            if(~istril(this.D))
                error('Diagonally Implicit Block Method must have a lower triangular D Matrix');
            end
            if(~isstril(this.C))
                error('Diagonally Implicit Block Method must have a strictly lower triangular C Matrix');
            end
        end
        
        function setNonZeroMatrixIndices(this)
        %SETNONZEROSTAGEINDICES set nonzero row entries of the block matrices A, B, C, D
        % = Parameters =================================================================================================
        % = Returns ====================================================================================================
        % ==============================================================================================================   
            this.non_zero_A_indices = this.nonzeroMatrixRows(this.A, true);
            this.non_zero_B_indices = this.nonzeroMatrixRows(this.B, true);
            this.non_zero_C_indices = this.nonzeroMatrixRows(this.C, false);
            this.non_zero_D_indices = this.nonzeroMatrixRows(this.D, false);            
        end

    end
    
end