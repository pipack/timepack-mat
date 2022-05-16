% =========================================================================
%
% Exponential polynomial block methods (EPBM) with Legendre nodes from [1]. 
%
%   options:
%       q     - number of nodes (default 3).
%       kappa - number of iterator applications (default 0).
%       alpha - extrapolation factor (default 2).
%
% [1] Buvoli, Tommaso. "Exponential polynomial block methods." SIAM Journal
% on Scientific Computing 43.3 (2021): A1692-A1722.
%
% =========================================================================

classdef epiPBMConst < IntegratorConst & ExponentialIntegratorConst
    
    properties
        graph_line_style = {};
        kappa = 0;
        alpha = 2;
    end
    
    properties(SetAccess = protected)
        name = '';
        description = '';
        order = [];
        starting_times   = []; 
        wts   = [];
    end
     
    properties(Access = private)
        q
        z
    end
    
    methods
        
        function this = epiPBMConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_field_value_pairs = {
                {'q', 3}
                {'kappa', 0}
                {'alpha', 2}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            this = this@IntegratorConst(options);
            this = this@ExponentialIntegratorConst(options);
            this.kappa = options.kappa;
            this.alpha = options.alpha;
            this.setQ(options.q);
        end
        
        function setQ(this, q)
            leg = legpts(q - 1, [0 2]);
            if(q > 3)
                this.order = q + 1;
            else
                this.order = q;
            end
            this.q               = q;
            this.name            = ['EPIPBM', num2str(q)];
            this.wts             = weights(0, leg, q - 2);
            this.starting_times  = 0;
            this.z               = [0, leg(:).'];           
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, t_in, y_in, problem)
                        
            ode_dim = size(y_in, 1);

            step_struct = struct( ...
                'r',    this.h / this.alpha, ....
                'b',    zeros(ode_dim, this.q), ...
                'rR',   zeros(ode_dim, this.q - 1), ...
                'rF_n', zeros(ode_dim, 1), ...
                'rJ_n', [] ...
            );
        
            y_in = repmat(y_in, [1, this.q]);
            
            % -- iterate for initial conditions ------------------------------------------------------------------------
            for i = 1 : this.q + this.kappa
                [~, y_in, step_struct] = this.pstep(0, t_in, y_in, step_struct, problem);
            end
                
        end
        
        function [t_out, y_out, step_struct, exit_flag] = pstep(this, alpha, t_in, y_in, step_struct, problem)
           
            % -- evaluate Remainder ---------------------------------------            
            step_struct = this.evalR(y_in, step_struct, problem);            
            % -- form vectors b -------------------------------------------
            step_struct.b(:, 2:end) = step_struct.rR * this.wts;
            step_struct.b(:, 2) = step_struct.b(:,2) + step_struct.rF_n;
            % -- advance solution --------------------------------------------------------------------------------------
            [phi_rR, exit_flag] = this.phi_evaluator.solve(this.z + alpha, step_struct.rJ_n, step_struct.b);
            y_out = y_in(:, 1) + phi_rR;
            t_out = t_in + step_struct.r * alpha;
                 
        end
        
        function [step_struct, r] = evalR(this, y_in, step_struct, problem) % Evaluate Remainder at Nodes 2, ..., q
            
            r   = this.h / this.alpha;
            y_n = y_in(:, 1);
            
            rF_n = r * problem.RHS(y_n);                                    % calculate right hand side
            rJ_n = r * problem.J(y_n);                                      % calculate Jacobian
            rR   = @(y) r * problem.RHS(y) - rF_n - rJ_n*(y - y_n);         % remainder function
           
            for j = 1 : this.q - 1 % evaluate at nodes z_2, ..., z_q
                step_struct.rR(:, j) = rR(y_in(:, j + 1));
            end
            
            % store r * F_n and r * J_n
            step_struct.rF_n = rF_n;
            step_struct.rJ_n = rJ_n; 

        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem)
            
            step_start_time = tic;
            
            [t_out, y_out, step_struct, exit_flag] = pstep(this, this.alpha, t_in, y_in, step_struct, problem);
            if(exit_flag == false)
                t_out = NaN; y_out = NaN; return;
            end
            
            for i = 1 : this.kappa
                [t_out, y_out, step_struct, exit_flag] = pstep(this, 0, t_out, y_out, step_struct, problem);
                if(exit_flag == false)
                    t_out = NaN; y_out = NaN; return;
                end
            end
                                    
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
        function [t_user, y_user] = userOutput(this, t_in, y_in, struct_in, t_out, y_out, struct_out, problem)
            t_user = t_out;
            y_user = y_out(:, 1);
        end
        
    end
    
end