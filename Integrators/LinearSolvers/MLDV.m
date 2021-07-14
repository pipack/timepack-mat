classdef MLDV < LinearSolver
    
    properties
        stats;
    end
    
	properties(SetAccess = protected)
        allows_matrix_free = false; % this class only allows a matrix A to be passed in as first argument of solve and solveBC
    end
    
    methods
        
        function this = MLDV(options)
            if(nargin == 0)
                options = struct();
            end
            options = this.DefaultOptions(options);
            this.stats = MLDVStats();
            this.stats.record = options.record_stats;
        end
        
        function [x, exit_flag, residual] = solve(this, A, b, varargin)            
            starting_time = tic;
            x = A\b;
            exit_flag = true;
            if(nargout > 2)
                residual = norm(A*x - b, 2);
            end            
            this.stats.recordSolve(toc(starting_time));
        end
        
        function [x, exit_flag, residual] = solveBC(this, A, b, c, varargin)
            % solves (I - c * A) x = b
            starting_time = tic;
            n = size(A,1);
            x = (speye(n) - c * A)\b;  % (I - c * A) = b 
            exit_flag = true;
            if(nargout > 2)
                residual = norm(A*x - b, 2);
            end
            this.stats.recordSolve(toc(starting_time));
        end
        
        function [x, exit_flag, residual] = solveBCKronA(this, A, b, C, varargin)
        % As (cell array of n matrices) 
        % b  (vector n * size(As,1))
        % C  (matrix nxn)
            
            % solves (I - kron(C,A)) x = b
            starting_time = tic;
            n_C = size(C,1);
            n_A = size(A,1);
            x = (speye(n_C * n_A) - kron(C, A)) \ b;
            exit_flag = true;
            if(nargout > 2)
                residual = norm(kron(C, A) * x - b, 2);
            end
            this.stats.recordSolve(toc(starting_time));
        end
        
        function [x, exit_flag, residual] = solveBCKronAs(this, As, b, C, varargin)
        % As (cell array of n matrices) 
        % b  (vector n * size(As,1))
        % C  (matrix nxn)
             
            % solves (I - kron(C,A)) x = b
            starting_time = tic;
            n_C = size(C,1);
            n_A = size(As{1},1);
            A_aug = speye(n_C * n_A) - kron(C, speye(n_A)) * blkdiag(As{:});
            x = A_aug \ b;
            exit_flag = true;
            if(nargout > 2)
                residual = norm(A_aug * x - b, 2);
            end
            this.stats.recordSolve(toc(starting_time));
        end
       
    end
    
    methods(Access = private)
        function options = DefaultOptions(this, options)
            if(~isfield(options, 'record_stats'))
                options.record_stats = false;
            end            
        end
    end
end