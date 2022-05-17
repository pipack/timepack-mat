classdef DI_RKConst < RKConst & ImplicitIntegratorConst
    
    properties
        graph_line_style = {};
    end
    
    properties(Abstract = true, SetAccess = protected)
        A % RK Matrix (Must be lower triangular)
        b % RK b vector
        c % RK c vector
    end
    
    properties(Abstract = true)
        eval_RHS % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
    end
    
    properties(SetAccess = protected)
        non_zero_stage_indices = {};
        non_zero_output_indices = [];
        nearest_stage_indices = {};
    end
    
    methods
        
        function this = DI_RKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@RKConst(options);
            this = this@ImplicitIntegratorConst(options);
            this.non_zero_stage_indices  = this.nonzeroStageIndices(this.A);
            this.non_zero_output_indices = find(this.b);
            this.nearest_stage_indices   = this.nearestStageIndices(this.c);
            % -- set up stat objects -----------------------------------------------------------------------------------
            num_stages = size(this.A, 1);
            this.nonlinear_solver.stats.reset(num_stages);
            this.nonlinear_solver.linear_solver.stats.reset(num_stages);
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct] = initStepStruct(this, ~, y_in, ~)
            
            num_stages = size(this.A, 1);
            ode_dim    = length(y_in);
            
            step_struct = struct(                   ...
                'AT', transpose(this.A),             ...
                'bT', transpose(this.b),             ...
                's',  num_stages,                   ...
                'Y',  zeros(ode_dim, num_stages),   ... % Stage Vector
                'F',  zeros(ode_dim, num_stages)    ... % Stage Derivative Vector
                );
            
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, ~)
            
            step_start_time = tic; % -- start step time clock ----------------------------------------------------------
            
            h = this.h; % stepsize
            num_stages = step_struct.s;
            
            for j = 1 : num_stages
                
                this.nonlinear_solver.stats.setIndex(j);
                this.nonlinear_solver.linear_solver.stats.setIndex(j);
                
                nzi = this.non_zero_stage_indices{j};
                nsi = this.nearest_stage_indices{j};
                
                % -- newton solve Y_j = b_j + c_j * f(Y_j) -------------------------------------------------------------
                b_j   = y_in + step_struct.F(:, nzi) * (h * step_struct.AT(nzi, j));
                c_j   = h * step_struct.AT(j, j);
                if(nsi == 0) % use y_in as initial guess
                    step_struct.Y(:,j) = this.nonlinear_solver.solveBC(problem, b_j, c_j, y_in);
                else % use Y_nsi as initial guess
                    step_struct.Y(:,j) = this.nonlinear_solver.solveBC(problem, b_j, c_j, step_struct.Y(:, nsi));
                end
                
                if(this.eval_RHS || c_j == 0)
                	start_rhs_time = tic;
                    step_struct.F(:, j) = problem.RHS(step_struct.Y(:,j));
                    this.rhs_stats.recordRHSEval(toc(start_rhs_time));
                else
                    step_struct.F(:, j) = (step_struct.Y(:,j) - b_j) / c_j;
                end
            end
            
            nzo = this.non_zero_output_indices;
            
            y_out = y_in + step_struct.F(:, nzo) * (h * step_struct.bT(nzo));
            t_out = t_in + h;
            
            this.step_stats.recordStep(toc(step_start_time)); % -- stop step time clock --------------------------------
            
        end
        
    end
    
end