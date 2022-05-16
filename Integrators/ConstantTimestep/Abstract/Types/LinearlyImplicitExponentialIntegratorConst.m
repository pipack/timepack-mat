classdef LinearlyImplicitExponentialIntegratorConst < IntegratorConst
    
    properties
        linear_solver
        phi_evaluator
        % mode = 1 -> f^{1} linearly implicit, f^{2} exponential
        % mode = 2 -> f^{1} exponential, f^{2} linearly implicit
        mode    
    end
    
    methods
        
        function this = LinearlyImplicitExponentialIntegratorConst(options)
            % -- set basic integrator properties -----------------------------------------------------------------------
            this@IntegratorConst(options);
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'linear_solver', MLDV()},  ...
                {'phi_evaluator', KIOPS()}  ...
                {'mode', 1}                 ...
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            this.phi_evaluator = options.phi_evaluator;
            this.linear_solver = options.linear_solver;
            this.mode = 1;
        end
           
        function J = RHSli(this, problem, u)
            switch(this.mode)
                case 1
                    J = problem.RHS(u, 1);
                case 2
                    J = problem.RHS(u, 2);
            end  
        end
        
        function J = RHSexp(this, problem, u)
            switch(this.mode)
                case 1
                    J = problem.RHS(u, 2);
                case 2
                    J = problem.RHS(u, 1);
            end
        end
       
        function J = Jli(this, problem, u)
            switch(this.mode)
                case 1
                    J = problem.J(u, 1);
                case 2
                    J = problem.J(u, 2);
            end  
        end
        
        function J = Jexp(this, problem, u)
            switch(this.mode)
                case 1
                    J = problem.J(u, 2);
                case 2
                    J = problem.J(u, 1);
            end
        end
        
    end
    
end