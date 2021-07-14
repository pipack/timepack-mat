% ======================================================================================================================
% Wraps MATLAB GMRES Linear Solver
% ======================================================================================================================

classdef GMRESMAT < LinearSolver
    
    properties
        tolerance;          % (double) tolerance of linear solver. iteration will exit if residual norm is below tolerance.
        max_iterations;     % (double) max number of iterations for building krylov basis
        force_real;         % (bool) If true, then matrix systems with imaginary values will be represented as real systems of twice the dimension.
        stats;              % (KrylovStats) stats object
        preconditioners;    % cell array of preconditioners
        restart;            % if non-empty, GMRES will restart arnoldi every restart iterations.
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
                {'restart', []}
                {'preconditioners', {[], []}}
                };
            options = setDefaultOptions(options, default_options_cell);
            this.max_iterations = options.max_iterations;
            this.tolerance = options.tolerance;
            this.force_real = options.force_real;
            this.stats = KrylovStats();
            this.stats.record = options.record_stats;
            this.preconditioners = options.preconditioners;
            this.restart = options.restart;
        end
        
        function [y, exit_flag, residual] = solve(this, A, b, x0)
            
            if(nargin == 4)
                x0 = zeros(size(b));
            end 
            
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
                        [y_aug, flag, residual, iter, resvec] = gmres(A_aug, b_aug, [], this.tolerance, min(2*dim, this.max_iterations), this.preconditioners{:}, x0_aug);
                        y = y_aug(1:dim) + 1i*y_aug(dim+1:2*dim);
                    else
                        [y, flag, residual, iter, resvec] = gmres(A, b, this.restart, this.tolerance, min(dim, this.max_iterations), this.preconditioners{:}, x0);
                    end
                else % -- function handle ------------------------------------------------------------------------------
                    b_aug  = [real(b); imag(b)];
                    x0_aug = [real(x0); imag(x0)];
                    [y_aug, flag, residual, iter, resvec] = gmres(@A_aug_handle, b_aug, [], this.tolerance, min(2*dim, this.max_iterations), this.preconditioners{:}, x0_aug);
                    y = y_aug(1:end/2) + 1i * y_aug(end/2+1:end);
                end
            else
                [y, flag, residual, iter, resvec] = gmres(A, b, this.restart, this.tolerance, min(dim, this.max_iterations), this.preconditioners{:}, x0);
            end
       
            if(flag == 0)
                exit_flag = true;
            else
                exit_flag = false;
            end
            this.stats.recordSolve(iter, resvec, toc(starting_time));
        end
        
        function [y, exit_flag, residual] = solveBC(this, A, b, c, x0)
        % A - matrix or function handle
        % c - constant
        % b - vector of dimension size(A,1)
        % x0 - optional initial guess
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
        
        function [y, exit_flag, residual] = solveBCKronA(this, A, b, C, x0)
        % A - matrix or function handle
        % C - square matrix
        % b - vector of dimension (size(A,1) * size(C,2))
        % x0 - optional initial guess
            
            if(nargin == 4)
                x0 = zeros(size(b));
            end
            
            C_n = size(C,1);
            A_n = length(b) / C_n; % deduce from b (in case A is a handle)
            
            function p = kronProdHandle(x) % function handle that returns (I - kron(M,A))                    
                p = x;
                for i = 1 : C_n
                    i_inds = (1 : A_n) + (i - 1) * A_n;
                    for j = 1 : C_n
                        j_inds = (1 : A_n) + (j - 1) * A_n;
                        p(i_inds) = p(i_inds) - C(i,j) * A(x(j_inds));
                    end
                end               
            end
            
            if(isnumeric(A))
                A_hat = speye(A_n * C_n) - kron(C, A);
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            else
                [y, exit_flag, residual] = solve(this, @kronProdHandle, b, x0);
            end  

        end
        
        function [y, exit_flag, residual] = solveBCKronAs(this, As, b, C, x0)
        % As - cell of nxn matrices or cell of function handles 
        % C - square matrix
        % b - vector of dimension (size(A,1) * size(C,2))
        % x0 - optional initial guess
            
            if(nargin == 4)
                x0 = zeros(size(b));
            end
            
            C_n = size(C,1);
            A_n = length(b) / C_n; % deduce from b (in case A is a handle)
            
            function p = kronProdHandle(x) % function handle that returns (I - kron(C, speye(n_A)) * blkdiag(As{:}));                    
                p = x;
                for i = 1 : C_n
                    i_inds = (1 : A_n) + (i - 1) * A_n;
                    for j = 1 : C_n
                        j_inds = (1 : A_n) + (j - 1) * A_n;
                        p(i_inds) = p(i_inds) - C(i,j) * As{j}(x(j_inds));
                    end
                end               
            end
            
            if(isnumeric(As{1}))
                A_hat = speye(A_n * C_n) - kron(C, speye(A_n)) * blkdiag(As{:});
                [y, exit_flag, residual] = solve(this, A_hat, b, x0);
            else
                [y, exit_flag, residual] = solve(this, @kronProdHandle, b, x0);
            end  

        end
        
    end
        
end