classdef IMEXConst < ImplicitIntegratorConst & LinearlyImplicitIntegratorConst
    %IMEXCONST Implicit-Explicit Integrator
     
    properties(Abstract)
        linearly_implicit % if set to true, then only linear solver will be used for implicit component (triggers automatically if linearize = true) 
        linearize % if set to true implicit component will be \frac{\partial J}{\partial x}
    end
    
    methods

        function this = IMEXConst(options)
            this = this@ImplicitIntegratorConst(options);
            this = this@LinearlyImplicitIntegratorConst(options);
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {{'linearly_implicit', false}, {'linearize', false}};
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            this.linearly_implicit = options.linearly_implicit;
            this.linearize = options.linearize;
        end
        
    end
    
    methods(Access = protected)
        
        function f = evalF(this, problem, y, y0, part)
            % problem - object that defines ODE
            % y - where y should be evaluated
            % y0 - local linearization point (only used if this.linearize == true)
            % part - value of 1 implies (implict) and value of 2 implies explicit 
            
            if(this.linearize)
                switch part
                    case 1
                        f = problem.J(y0) * y; % full jacobian
                    case 2
                        f = problem.RHS(y) - problem.J(y0) * y; % full RHS - full jacobian
                end
            else
                f = problem.RHS(y, part);
            end
        end
        
        function j = evalJ(this, problem, y0, part) % component jacobian
            if(this.linearize)
                switch part
                    case 1
                        j = problem.J(y0); % full jacobian
                    case 2
                        dim = length(y0);
                        j   = sparse([],[],[], dim, dim);
                end
            else
                j = problem.J(y0, part);
            end
        end
        
        function [x, clean_exit] = solveBC(this, problem, b, c, y_guess, y0)
            % SOLVEBC solves linear or nonlinear problem associated
            %
            % If this.linearly_implicit = true, solves the linear system
            %
            %   (I - c * problem.J(y0)) * x = b
            %
            % If this.linearly_implicit = false, solve the nonlinear system
            %
            %   (x - c * problem.RHS(x)) = b
            %
            % -- Paramters --------------------------------
            %  problem - object that defines ode
            %  b (vector)
            %  c (constant)
            %  y_guess (vector) - guess for linear or nonlinear solver
            %  y0 - location of linearization (only used if this.linearly_implicit = true)
            
            part = 1; % treat first part implicitly
            if(this.linearly_implicit || this.linearize)
               [x, clean_exit] = this.linear_solver.solveBC(this.evalJ(problem, y0, part), b, c, y_guess);
            else
               [x, clean_exit] = this.nonlinear_solver.solveBC(problem, b, c, y_guess, part);
            end
            
        end
        
        function [x, clean_exit] = solveBCKronF(this, problem, b, C, y_guess, y0)
            % SOLVEBC solves linear or nonlinear problem associated
            %
            % If this.linearly_implicit = true, solves the linear system
            %
            %   (I - kron(C,problem.J(y0)) * x = b
            %
            % If this.linearly_implicit = false, solve the nonlinear system
            %
            %   x - kron(C, F(x)) = b 
            %
            % -- Paramters --------------------------------
            %  problem - object that defines ode
            %  b (vector)
            %  c (constant)
            %  y_guess (vector) - guess for linear or nonlinear solver
            %  y0 - location of linearization (only used if this.linearly_implicit = true)
            
            part = 1; % treat first part implicitly
            if(this.linearly_implicit || this.linearize)
               [x, clean_exit] = this.linear_solver.solveBCKronA(this.evalJ(problem, y0, part), b, C, y_guess);
            else
               [x, clean_exit] = this.nonlinear_solver.solveBCKronF(problem, b, C, y_guess, part);
            end
            
        end

    end
    
end