% =========================================================================
%
% 2nd-order Exponential Euler Method
%
% =========================================================================

classdef epi2Const < IntegratorConst & ExponentialIntegratorConst
    
    properties
        graph_line_style = {};
        matrix_free = false;
    end
    
    properties(SetAccess = protected)
        name = '';
        description = '';
        order = [];
        starting_times = [0];
    end
    
    methods
        function this = epi2Const(options)
            if(nargin == 0)
                options = struct();
            end
            default_field_value_pairs = {};
            options = setDefaultOptions(options, default_field_value_pairs);
            this = this@IntegratorConst(options);
            this = this@ExponentialIntegratorConst(options);
            
            this.order = 2;
            this.name  = 'EPI2';
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
            step_struct = struct();
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, final_step)
            
            step_start_time = tic;
            
            % -- store parameters into local variables --------------------
            h = this.h;
            % -- evaluate right-hand-side ---------------------------------
            hF_n = h * problem.RHS(y_in);  % calculate right hand side
            if(this.matrix_free)
                hJ_n = @(x) h * problem.Jx(y_in, x);  % function handle for product
            else
                hJ_n = h * problem.J(y_in);  % calculate full Jacobian
            end
            % -- advance solution -----------------------------------------
            u = [zeros(length(hF_n),1), hF_n];
            y_out = y_in + this.phi_evaluator.solve(1, hJ_n, u);
            t_out = t_in + h;
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
    end
    
end