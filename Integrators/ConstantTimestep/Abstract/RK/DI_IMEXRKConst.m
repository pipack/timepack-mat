% ======================================================================================================================
%
%  Diagonally Implicit - Explicit RK method for solving
%
%       y' = fi(y) + fe(y)
%
%   The first term fe(y) will be treated implicitly while the second term fe(y) will be treated explicitly.
%
% ======================================================================================================================

classdef DI_IMEXRKConst < RKConst & IMEXConst
    
    properties
        ilu_precondition = false % if true and linear_solver is of type GMRESMAT, then an ILU preconditioner for (I - h * D(1,1) * L ) will be precomputed when stepsize is set
    end
    
    properties(Abstract = true, SetAccess = protected)
        % implicit component
        Ai % RK Matrix (must be lower triangular)
        bi % RK b vector
        % explicit component
        Ae % RK Matrix (explicit component, must be strictly lower triangular)       
        be % RK b vector
        c % RK c vector
    end
    
    properties(Abstract = true)
        eval_RHS % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
    end
    
    properties(SetAccess = protected)
        non_zero_stage_indices_Ai = {};
        non_zero_stage_indices_Ae = {};
        non_zero_output_indices_bi = [];
        non_zero_output_indices_be = [];
        nearest_stage_indices = {};
    end
    
    methods
        
        function this = DI_IMEXRKConst(options)
            if(nargin == 0)
                options = struct();
            end
            options = setDefaultOptions(options, {{'ilu_precondition', false}});
            this = this@RKConst(options);
            this = this@IMEXConst(options);
            this.non_zero_stage_indices_Ai  = this.nonzeroStageIndices(this.Ai);
            this.non_zero_stage_indices_Ae  = this.nonzeroStageIndices(this.Ae);
            this.non_zero_output_indices_bi = find(this.bi);
            this.non_zero_output_indices_be = find(this.be);
            this.nearest_stage_indices   = this.nearestStageIndices(this.c);
            % -- set up stat objects -----------------------------------------------------------------------------------
            num_stages = size(this.Ai, 1);
            this.linear_solver.stats.reset(num_stages);
            % precoditioner
            this.ilu_precondition = options.ilu_precondition;
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct] = initStepStruct(this, ~, y_in, ~)
            
            num_stages = size(this.Ai, 1);
            ode_dim    = length(y_in);
            
            step_struct = struct(                   ...
                'AiT', transpose(this.Ai),          ...
                'biT', transpose(this.bi),          ...
                'AeT', transpose(this.Ae),          ...
                'beT', transpose(this.be),          ...
                's',   num_stages,                   ...
                'Y',   zeros(ode_dim, num_stages),   ... % Stage Vector
                'Fi',  zeros(ode_dim, num_stages),    ... % Stage Derivative Vector
                'Fe',  zeros(ode_dim, num_stages)    ... % Stage Derivative Vector
                );
            
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, ~)
            
            step_start_time = tic; % -- start step time clock ----------------------------------------------------------
            
            h = this.h; % stepsize
            num_stages = step_struct.s;
            
            for j = 1 : num_stages
                
                % -- init stats ----------------------------------------------------------------------------------------
                if(this.linearly_implicit)
                    this.linear_solver.stats.setIndex(j);
                else
                    this.nonlinear_solver.stats.setIndex(j);
                    this.nonlinear_solver.linear_solver.stats.setIndex(j);
                end
                
                nzi_Ai = this.non_zero_stage_indices_Ai{j};
                nzi_Ae = this.non_zero_stage_indices_Ae{j};
                nsi = this.nearest_stage_indices{j};
                
                % -- solve Y_j = b_j + c_j * f(Y_j) --------------------------------------------------------------------
                b_j   = y_in + step_struct.Fi(:, nzi_Ai) * (h * step_struct.AiT(nzi_Ai, j)) + step_struct.Fe(:, nzi_Ae) * (h * step_struct.AeT(nzi_Ae, j));
                c_j   = h * step_struct.AiT(j, j);
                
                if(c_j ~= 0)                
                    if(nsi == 0) % use y_in as initial guess
                        step_struct.Y(:,j) = this.solveBC(problem, b_j, c_j, y_in, y_in);
                    else % use Y_nsi as initial guess
                        step_struct.Y(:,j) = this.solveBC(problem, b_j, c_j, step_struct.Y(:, nsi), y_in);
                    end
                else
                    step_struct.Y(:,j) = b_j;
                end
                
                if(this.eval_RHS || c_j == 0)
                	start_rhs_time = tic;
                    step_struct.Fi(:, j) = this.evalF(problem, step_struct.Y(:,j), y_in, 1);
                    this.rhs_stats.recordRHSEval(toc(start_rhs_time));
                else
                    step_struct.Fi(:, j) = (step_struct.Y(:,j) - b_j) / c_j;
                end
                
                start_rhs_time = tic;
                step_struct.Fe(:, j) = this.evalF(problem, step_struct.Y(:,j), y_in, 2);
                this.rhs_stats.recordRHSEval(toc(start_rhs_time));
                
            end
            
            nzo_bi = this.non_zero_output_indices_bi;
            nzo_be = this.non_zero_output_indices_be;
            
            y_out = y_in + step_struct.Fi(:, nzo_bi) * (h * step_struct.biT(nzo_bi)) + step_struct.Fe(:, nzo_be) * (h * step_struct.beT(nzo_be));
            t_out = t_in + h;
            
            this.step_stats.recordStep(toc(step_start_time)); % -- stop step time clock --------------------------------
            
        end
        
        function setStepsize(this, problem)
            setStepsize@IntegratorConst(this, problem)
            
            nz_diag_coeff = setdiff(unique(diag(this.Ai)), 0);
            if(this.ilu_precondition && ~this.linearize && isa(this.linear_solver, 'GMRESMAT') && length(nz_diag_coeff) == 1)
                [L, U] = ilu(speye(problem.dimension) - this.h * nz_diag_coeff(1) * this.evalJ(problem, problem.initial_condition, 1));
                this.linear_solver.preconditioners = {L, U};
            end            
        end
        
    end
    
end