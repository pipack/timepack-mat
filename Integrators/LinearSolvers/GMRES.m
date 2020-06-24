% ======================================================================================================================
% Classical GMRES Linear Solver
% ======================================================================================================================

classdef GMRES < LinearSolver
    
    properties
        tolerance;      % (double) tolerance of linear solver. iteration will exit if residual norm is below tolerance.
        max_iterations  % (double) max number of iterations for building krylov basis
        tol_norm        % (handle) norm used to determine residual exit tolerance 
        force_real;     % (bool) If true, then matrix systems with imaginary values will be represented as real systems of twice the dimension.
        stats;          % (KrylovStats) stats object
    end
    
	properties(SetAccess = protected)
        allows_matrix_free = true; % this class allows both a matrix A or a function handles that computes A*x to be passed in as first argument of solve and solveBC 
    end
    
    properties (Access = private)
        m;
        vec_length;
        H;
        H2;
        V;
    end
    
    methods
        
        function this = GMRES(options)
            if(nargin == 0)
                options = struct();
            end
            default_options_cell = {
                {'tolerance', 1e-10}
                {'max_iterations', 100}
                {'tol_norm', @(x) norm(x, 2)}
                {'record_stats', false}
                {'force_real', false}
                };
            options = setDefaultOptions(options, default_options_cell);
            this.setMaxIter(options.max_iterations);
            this.tolerance = options.tolerance;
            this.tol_norm = options.tol_norm;
            this.force_real = options.force_real;
            this.stats = KrylovStats();
            this.stats.record = options.record_stats;
        end
        
        function setMaxIter(this, m)
            this.m  = m;
            this.H  = zeros(m, m);
        end
        
        function [y, exit_flag, residual] = solve(this, A, b, x0)
            
            dim = length(x0);
            function y = A_aug_handle(x)
                y = A(x(1:dim) + 1i * x(dim+1:2*dim));
                y = [real(y); imag(y)];
            end
            
            if(this.force_real)
                if(isnumeric(A)) % explicit jacobian matrix ------------------------------------------------------------
                    real_inputs = (isreal(A) && isreal(b) && isreal(x0));
                    if(~real_inputs)
                        A_r = real(A);
                        A_i = imag(A);
                        A_aug  = [A_r -A_i; A_i A_r];
                        b_aug  = [real(b); imag(b)];
                        x0_aug = [real(x0); imag(x0)];
                        [y_aug, exit_flag, residual] = this.gmres(A_aug, b_aug, x0_aug);
                        y = y_aug(1:dim) + 1i*y_aug(dim+1:2*dim);
                    else
                        [y, exit_flag, residual] = this.gmres(A, b, x0);
                    end
                else % -- function handle ------------------------------------------------------------------------------
                    b_aug  = [real(b); imag(b)];
                    x0_aug = [real(x0); imag(x0)];
                    [y_aug, exit_flag, residual] = this.gmres(@A_aug_handle, b_aug, x0_aug);
                    y = y_aug(1:end/2) + 1i * y_aug(end/2+1:end);
                end
            else
                [y, exit_flag, residual] = this.gmres(A, b, x0);
            end  
        end
        
        
        function [y, exit_flag, residual] = solveBC(this, A, b, c, x0) % solves x = b + c * A * x           
            if(isnumeric(A))
                A_hat = speye(size(A)) - c * A;
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            else
                A_hat = @(x) x - c * A(x);
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            end
            
        end
        
    end
    
    methods(Access = private)
        
        function [x, exit_flag, residual, j, Vout, Hout] = gmres(this, A, b, x0)
            
            verifyDimension(this, length(x0));
            if(isnumeric(A)) % assume A is matrix
                A_is_numeric = true;
            else % Assume A is function handle
                A_is_numeric = false;
            end
            
            % -- initialize stats --
            residuals = zeros(this.m + 1, 1);
            seconds   = zeros(this.m, 1);
            exit_flag = true;
            
            % Pre-compute the norm of the r0 vector.  If it's too small, exit.
            if(A_is_numeric)
                r0 = b - A*x0;
            else
                r0 = b - A(x0);
            end
            
            beta = norm(r0);
            if beta < eps  % if beta below machine precision
                % Essentially beta is a zero vector, there's nothing to do.
                Vout = r0;
                Hout = 1;
                x = x0;
                residuals(1) = beta; 
                residual = beta;               
                this.stats.recordSolve(0, residuals(1), 0); 
                return;
            end
            
            % Arnoldi iteration.
            maxBasis = this.m;
            this.V(:,1) = r0 / beta;
            for j = 1:maxBasis
                iteration_start_time = tic;
                
                if(A_is_numeric)
                    w = A * this.V(:,j);
                else
                    w = A(this.V(:,j));
                end
                for i = 1:j
                    this.H(i,j) = this.V(:,i)' * w;
                    w = w - this.H(i,j) * this.V(:,i);
                end
                this.H(j+1,j) = norm(w);
                this.V(:,j+1) = w / this.H(j+1,j);                
                y = this.H(1:j+1,1:j) \ (beta*eye(j+1,1));
                
                residuals(j+1) = this.tol_norm(this.H(1:j+1,1:j)*y - beta*eye(j+1,1));
                if (residuals(j+1) < this.tolerance) % tolerance achieved. exit with exit_flag = true
                    Hout = this.H(1:j, :);
                    Vout = this.V(:, 1:j);
                    x = x0 + this.V(:, 1:j) * y;
                    residual = residuals(j+1);                    
                    seconds(j) = toc(iteration_start_time);
                    this.stats.recordSolve(j, residuals(1:j+1), seconds(1:j));
                    return;
                end
                
                seconds(j) = toc(iteration_start_time);
            end
            
            warning('GMRES exceeded max iterations; final residual was: %e.', residuals(j));
            exit_flag = false;
            % Truncate off extra final iteration.
            Hout = this.H(1:j, :);
            Vout = this.V(:, 1:j);
            x = x0 + this.V(:, 1:j) * y;
            residual = residuals(j);
            
            % -- STORE KRYLOV STATS --
            this.stats.recordSolve(j, residuals(1:j+1), seconds(1:j));            
            
        end        
        
        function verifyDimension(this, vec_length)
            if(isempty(this.vec_length))
                reset_flag = true;
            elseif(this.vec_length ~= vec_length)
                reset_flag = true;
            else
                reset_flag = false;
            end
            if(reset_flag)
                this.vec_length = vec_length;
                this.V = zeros(vec_length, this.m);
            end
        end       
        
    end
    
end