% ======================================================================================================================
% Newton Nonlinear Solver
% ======================================================================================================================

classdef Newton < NonlinearSolver
    
    properties
        max_iterations              % (double) max number of iterations before newton exits
        max_residual_tolerance      % (double) newton's method will exit with clean_exit_flag = false if residual is greater than this quantity
        residual_tolerance          % (double) newton's method will exit with clean_exit_flag = true if residual is below this quantity
        delta_tolerance             % (double) newton's method will exit with clean_exit_flag = true if change in residual compared to previous step is below this quantity
        tol_norm                    % (handle)
        linear_solver               % (LinearSolver) linear solver to apply at each iteration
        matrix_free                 % (bool) if true uses matrix-vector products problem.Jx otherwise uses matrix from problem.J
        stats                       % (NewtonStats) stats object
    end
    
    methods
        
        function this = Newton(options)
            if(nargin == 0)
                options = struct();
            end
            default_options = {
                {'max_iterations', 100}
                {'residual_tolerance', 1e-12}
                {'max_residual_tolerance', Inf}
                {'delta_tolerance', 1e-12}
                {'tol_norm', @(x) norm(x, Inf)}
                {'linear_solver', MLDV()}
                {'matrix_free', false}
                {'record_stats', false}
                };
            options = setDefaultOptions(options, default_options);
            this.stats = NewtonStats();
            this.max_iterations = options.max_iterations;
            this.residual_tolerance = options.residual_tolerance;
            this.max_residual_tolerance = options.max_residual_tolerance;
            this.delta_tolerance = options.delta_tolerance;
            this.tol_norm = options.tol_norm;
            this.linear_solver = options.linear_solver;
            this.matrix_free = options.matrix_free;
            this.stats.record = options.record_stats;
        end
        
        function set.matrix_free(this, value)
            if(value == true && ~this.linear_solver.allows_matrix_free) % do not allow matrix_free flag with MLDV
                this.matrix_free = false;
                warning('matrix_free newton not allowed with linear solver of type %s. Flag was set to false', class(this.linear_solver));
            else
                this.matrix_free = value;
            end
        end
        
        function [x_k, clean_exit_flag, final_residual] = solve(this, problem, x0, part)
            % SOLVEBC solves system 0 = F(x)
            % PARAMETERS
            %   problem (class or struct)   : must have the parameters:
            %                                   varArgs (Cell)   - additional arguments for F and J
            %                                   F       (Handle) - function for F(x)
            %                               and one of the following
            %                                   Jxy (Handle) - function @(x,y) returning product of jacobian at x with vector y
            %                                   J (Handle) - function @(x) returning matrix representing the jacobian of F at x
            %   b       (vector)            : constant b in nonlinear system
            %   c       (double)            : constant c in nonlinear system
            %   x0      (vector)            : initial guess
            
            if(nargin == 3)
                part_args = {};
            else
                part_args = {part};
            end
            
            % -- define function ---------------------------------------------------------------------------------------
            function Gx = G(x)
                Gx = problem.RHS(x, part_args{:});
            end
            
            % -- define Jacobian ---------------------------------------------------------------------------------------
            function GPx = Gp(x) % matrix formulation (returns matrix)
                GPx  = problem.J(x, part_args{:});
            end

            function GPxH = GpMF(x) % matrix free formulation (returns handle)
                GPxH = @(y) problem.Jx(x, y, part_args{:});
            end
            
            % -- Newton solve ------------------------------------------------------------------------------------------
            if( this.matrix_free )
                [x_k, clean_exit_flag, final_residual] = this.newtonSolver(@G, @GpMF, x0);
            else
                [x_k, clean_exit_flag, final_residual] = this.newtonSolver(@G, @Gp, x0);
            end
            
        end
        
        function [x_k, clean_exit_flag, final_residual] = solveBC(this, problem, b, c, x0, part)
        % SOLVEBC solves system x = b + c * F(x)
        % PARAMETERS
        %   problem (class or struct)   : must have the parameters:
        %                                   varArgs (Cell)   - additional arguments for F and J
        %                                   F       (Handle) - function for F(x)
        %                               and one of the following
        %                                   Jxy (Handle) - function @(x,y) returning product of jacobian at x with vector y
        %                                   J (Handle) - function @(x) returning matrix representing the jacobian of F at x
        %   b       (vector)            : constant b in nonlinear system
        %   c       (double)            : constant c in nonlinear system
        %   x0      (vector)            : initial guess
        %   part    (integer)           : options, specify if only F is only a part of the RHS
            if(nargin == 5)
                part_args = {};
            else
                part_args = {part};
            end
            
            % -- define function ---------------------------------------------------------------------------------------
            function Gx = G(x)
                Gx = b + c * problem.RHS(x, part_args{:}) - x;
            end
            
            % -- define Jacobian ---------------------------------------------------------------------------------------
            spI = [];
            if( ~ this.matrix_free )
                spI = speye(length(x0));
            end
            
            function GPx = Gp(x) % matrix formulation (returns matrix)
                GPx  = c * problem.J(x, part_args{:}) - spI;
            end

            function GPxH = GpMF(x) % matrix free formulation (returns handle)
                GPxH = @(y) c * problem.Jx(x, y, part_args{:}) - y;
            end
            
            % -- Newton solve ------------------------------------------------------------------------------------------
            if( this.matrix_free )
                [x_k, clean_exit_flag, final_residual] = this.newtonSolver(@G, @GpMF, x0);
            else
                [x_k, clean_exit_flag, final_residual] = this.newtonSolver(@G, @Gp, x0);
            end

        end
          
    end

    methods ( Access = private )
        
        function [x_k, clean_exit_flag, final_residual] = newtonSolver(this, G, GP, x0)
            % SOLVEBC solves system x = b + c * F(x)
            % PARAMETERS
            %   problem (class or struct)   : must have the parameters:
            %                                   varArgs (Cell)   - additional arguments for F and J
            %                                   F       (Handle) - function for F(x)
            %                               and one of the following
            %                                   Jxy (Handle) - function @(x,y) returning product of jacobian at x with vector y
            %                                   J (Handle) - function @(x) returning matrix representing the jacobian of F at x
            %   b       (vector)            : constant b in nonlinear system
            %   c       (double)            : constant c in nonlinear system
            %   x0      (vector)            : initial guess
            
            % -- Newton Iteration --------------------------------------------------------------------------------------
            iterations = 0;
            x_k        = x0;
            residuals  = zeros(this.max_iterations + 1, 1);
            deltas     = zeros(this.max_iterations, 1);
            seconds    = zeros(this.max_iterations, 1);
            clean_exit_flag = false;
            
            while(true) % DOWHILE loop
                iteration_start_time = tic;
                % -- Evaluate Nonlinear Function -----------------------------------------------------------------------
                G_k = G(x_k);
                % -- Calculate residual and test exit conditions -------------------------------------------------------
                residual_index = iterations + 1;
                current_residual = this.tol_norm(G_k);
                residuals(residual_index) = current_residual;
                if (current_residual < this.residual_tolerance) % clean residual exit condition
                    clean_exit_flag = true;
                    break;
                end
                if(iterations >= this.max_iterations || current_residual >= this.max_residual_tolerance || isnan(current_residual)) % bad exit conditions
                    break
                end
                % -- Construct Jacobian and Solve ----------------------------------------------------------------------
                GP_k = GP(x_k);
                
                if(isnumeric(GP_k))
                    rhs = GP_k * x_k - G_k;
                else
                    rhs = GP_k(x_k) - G_k;
                end
                
                x_km1 = x_k;
                x_k   = this.linear_solver.solve(GP_k, rhs, x_k);
                % -- increment iteration count and test exit conditions ------------------------------------------------
                iterations = iterations + 1;
                deltas(iterations) = this.tol_norm(x_k - x_km1);
                seconds(iterations) = toc(iteration_start_time);
                if(deltas(iterations) < this.delta_tolerance) % clean delta exit condition
                    clean_exit_flag = true;
                    break;
                end
                
            end
            final_residual = residuals(residual_index);
            this.stats.addSolve(iterations, residuals(1:iterations+1), deltas(1:iterations), seconds(1:iterations))
        end
         
    end

end