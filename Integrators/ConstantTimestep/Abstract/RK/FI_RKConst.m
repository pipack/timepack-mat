classdef FI_RKConst < RKConst & ImplicitIntegratorConst
    
    properties
        graph_line_style = {};
        
    end
    
    properties(Abstract = true, SetAccess = protected)
        A % RK Matrix (Must be lower triangular)
        b % RK b vector
        c % RK c vector
        A_extrapolate % initial guess extrapolation coefficients (matrix of dimension (s+2)xs input, outputs, stages)
    end
    
    properties(Abstract = true)
        eval_RHS % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
        extrapolate_initial_guess % initial guess extrapolation coefficients (matrix of dimension (s+2)xs input, outputs, stages)
    end
    
    properties(SetAccess = protected)
        non_zero_output_indices = [];
    end
    
    methods
        
        function this = FI_RKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@RKConst(options);
            this = this@ImplicitIntegratorConst(options);
            this.non_zero_output_indices = find(this.b);
            % -- set up stat objects -----------------------------------------------------------------------------------
            num_stages = size(this.A, 1);
            this.nonlinear_solver.stats.reset(num_stages);
            this.nonlinear_solver.linear_solver.stats.reset(num_stages);
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, ~, y_in, ~)
            
            num_stages = size(this.A, 1);
            ode_dim    = length(y_in);
            
            step_struct = struct(                   ...
                'ATinv', transpose(inv(this.A)),    ...
                'bT', transpose(this.b),            ...
                's',  num_stages,                   ...
                'Y',  zeros(ode_dim, num_stages),   ... % Stage Vector
                'F',  zeros(ode_dim, num_stages),   ... % Stage Derivative Vector
                'A_extrapT', transpose(this.A_extrapolate), ...
                'Y_guess', repmat(y_in, num_stages, 1) ...
            );
            
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, ~)
            
            step_start_time = tic; % -- start step time clock ----------------------------------------------------------
            h = this.h;
            
            % --> solve fully implicit nonlinear system
            Yn = repmat(y_in, step_struct.s, 1);
            if(this.extrapolate_initial_guess)
                Y_guess = step_struct.Y_guess;
            else
                Y_guess = Yn;
            end
            
            step_struct.Y = reshape( ...
                this.nonlinear_solver.solveBCKronF(problem, Yn, h * this.A, Y_guess), ... % solve x = bv + kron(h * A, F)
                size(step_struct.Y) ...
            );
            
            % --> evaluate RHS
            nzo = this.non_zero_output_indices;
            if(this.eval_RHS)
                for i = nzo(:)'
                    step_struct.F(:,i) = problem.RHS(step_struct.Y(:,i));
                end                
            else
                step_struct.F = ((step_struct.Y - reshape(Yn, size(step_struct.Y))) * step_struct.ATinv) / h;
            end
            
            % --> advance step            
            y_out = y_in + step_struct.F(:, nzo) * (h * step_struct.bT(nzo));
            t_out = t_in + h;
            
            if(this.extrapolate_initial_guess)
                A_extrapT = step_struct.A_extrapT;
                step_struct.Y_guess = reshape( ...
                    y_in * A_extrapT(1, :) ...
                    + step_struct.Y * A_extrapT(2:step_struct.s+1, :) ...
                    + y_out * A_extrapT(step_struct.s+2, :), ...
                    [numel(step_struct.Y),1] ...
                );
            end
                        
            this.step_stats.recordStep(toc(step_start_time)); % -- stop step time clock --------------------------------
            
        end
        
    end
    
end