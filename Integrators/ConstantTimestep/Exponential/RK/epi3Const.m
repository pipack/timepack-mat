classdef epi3Const < IntegratorConst & ExponentialIntegratorConst
    
    properties
        graph_line_style = {};
        matrix_free = false;
    end
    
    properties(SetAccess = protected)
        name = '';
        description = '';
        order = [];
        starting_times = [0 1];
        weights = [];
    end
    
    methods
        function this = epi3Const(options)
            if(nargin == 0)
                options = struct();
            end
            default_field_value_pairs = {};
            options = setDefaultOptions(options, default_field_value_pairs);
            this = this@IntegratorConst(options);
            this = this@ExponentialIntegratorConst(options);
            
            this.order = 3;
            this.name  = 'EPI3';
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct] = initStepStruct(this, t_in, y_in, problem)
            step_struct = struct(...
                'r', zeros(size(y_in,1), 2)    ...
                );
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem)
            
            step_start_time = tic;
            mtrx_free = this.matrix_free;
            
            % -- store parameters into local variables --------------------
            h = this.h;
            % -- evaluate remainder terms ---------------------------------
            y_n  = y_in(:,2);
            hF_n = h * problem.RHS(y_n);  % calculate right hand side
            if(mtrx_free)
                hJ_n = @(x) h * problem.Jx(y_n, x);  % function handle for product
            else
                hJ_n = h * problem.J(y_n);  % calculate full Jacobian
            end
            
            step_struct.r(:) = 0;
            if(mtrx_free)
                step_struct.r(:, 1) = h * problem.RHS(y_in(:, 1)) - hF_n - hJ_n(y_in(:, 1) - y_n);
            else
                step_struct.r(:, 1) = h * problem.RHS(y_in(:, 1)) - hF_n - hJ_n*(y_in(:, 1) - y_n);
            end
            % -- advance solution -----------------------------------------
            u = [zeros(size(step_struct.r,1),1), hF_n, 2/3 * step_struct.r(:,1)];
            y_out = y_n + this.phi_evaluator.solve([0 1], hJ_n, u);
            t_out = t_in + h;
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
        function [t_user, y_user] = userOutput(this, t_in, y_in, struct_in, t_out, y_out, struct_out, problem)
            t_user = t_out(1);
            y_user = y_out(:, 1);
        end
        
    end
    
end