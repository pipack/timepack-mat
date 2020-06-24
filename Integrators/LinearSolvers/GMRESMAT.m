% ======================================================================================================================
% Wraps MATLAB GMRES Linear Solver
% ======================================================================================================================

classdef GMRESMAT < LinearSolver
    
    properties
        tolerance;          % (double) tolerance of linear solver. iteration will exit if residual norm is below tolerance.
        max_iterations;     % (double) max number of iterations for building krylov basis
        force_real;         % (bool) If true, then matrix systems with imaginary values will be represented as real systems of twice the dimension.
        stats;              % (KrylovStats) stats object
    end
    
    properties(SetAccess = protected)
        allows_matrix_free = true; % this class allows both a matrix A or a function handles that computes A*x to be passed in as first argument of solve and solveBC 
    end
    
    methods
        
        function this = GMRESMAT(options)
            if(nargin == 0)
                options = struct();
            end
            default_options_cell = {
                {'tolerance', 1e-10}
                {'max_iterations', 100}
                {'record_stats', false}
                {'force_real', false}
                };
            options = setDefaultOptions(options, default_options_cell);
            this.max_iterations = options.max_iterations;
            this.tolerance = options.tolerance;
            this.force_real = options.force_real;
            this.stats = KrylovStats();
            this.stats.record = options.record_stats;
        end
        
        function [y, exit_flag, residual] = solve(this, A, b, x0)
            
            starting_time = tic;
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
                        [y_aug, flag, residual, iter, resvec] = gmres(A_aug, b_aug, [], this.tolerance, min(2*dim, this.max_iterations), [], [], x0_aug);
                        y = y_aug(1:dim) + 1i*y_aug(dim+1:2*dim);
                    else
                        [y, flag, residual, iter, resvec] = gmres(A, b, [], this.tolerance, min(dim, this.max_iterations), [], [], x0);
                    end
                else % -- function handle ------------------------------------------------------------------------------
                    b_aug  = [real(b); imag(b)];
                    x0_aug = [real(x0); imag(x0)];
                    [y_aug, flag, residual, iter, resvec] = gmres(@A_aug_handle, b_aug, [], this.tolerance, min(2*dim, this.max_iterations), [], [], x0_aug);
                    y = y_aug(1:end/2) + 1i * y_aug(end/2+1:end);
                end
            else
                [y, flag, residual, iter, resvec] = gmres(A, b, [], this.tolerance, min(dim, this.max_iterations), [], [], x0);
            end
       
            if(flag == 0)
                exit_flag = true;
            else
                exit_flag = false;
            end
            this.stats.recordSolve(iter, resvec, toc(starting_time));
        end
        
        
        function [y, exit_flag, residual] = solveBC(this, A, b, c, x0)
            if(nargin == 4)
                x0 = zeros(size(b));
            end
            if(isnumeric(A))
                A_hat = speye(size(A)) - c * A;
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            else
                A_hat = @(x) x - c * A(x);
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            end
            
        end
        
    end
        
end