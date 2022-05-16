
% ======================================================================================================================
%
%  Abstract Block Method:
%
%       y^{(n+1)} = A * y^{(n)} + r * B * y^{(n+1)} + C * y^{(n+1)} + r * D * f^{[n+1]}
%
% ======================================================================================================================

classdef BLKConst < IntegratorConst
    
    properties(Abstract = true, SetAccess = protected)
        % -- congugate properties --------------------------------------------------------------------------------------
        conjugate_inputs       % integer vector. if ith position is zero, then jth output should be computed. if jth position is non-zero then jth output has conjugate whose index is conjugate_inputs(j)  
        conjugate_outputs      % integer vector. if ith position is zero, then jth output should be computed. if jth position is non-zero then jth output has conjugate whose index is conjugate_outputs(j)
        real_valued_outputs     % bool vector. if ith position is true, then ith output can be computed using only real arithmetic after clipping imaginary parts
    end
    
	properties(Abstract = true)
    	eval_RHS                    % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
        parallel_initial_guess      % if false, previously computed outputs may be used as initial conditions. If true, only initial conditions will be used as guess 
        extrapolate_initial_guess   % if true coefficients A_extrapolate and B_extrapolate will be used to form initial guess for any nonlinear systems
    end
    
    properties(SetAccess = protected)
        starting_times = 0;
    end
    
    methods
    
        function this = BLKConst(options)  % -- Constructor ------------------------------------------------------------
            default_field_value_pairs = { ...
                {'extrapolate_initial_guess', true} ...
                {'eval_RHS', false} ...
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            this@IntegratorConst(options); % set basic integrator properties
            this.eval_RHS = options.eval_RHS;
            this.extrapolate_initial_guess = options.extrapolate_initial_guess;            
        end
                
    end
    
    methods (Access = protected) % -- generic helper functions ---------------------------------------------------------
                
        function indices = nonzeroMatrixRows(~, A, include_future_outputs)
        %NONZEROSTAGEINDICES compute nonzero row entries of the matrix A
        % = Parameters =================================================================================================
        %   1. A (Matrix) : matrix to get nonzero rows
        % = Returns ====================================================================================================
        %   1. indices (cell) : cell array where indices{i} contains indices of A that are non-zero
        % ==============================================================================================================
            num_outputs = size(A, 1);
            indices     = cell(num_outputs, 1);
            if(include_future_outputs)            
                for i = 1 : num_outputs
                    indices{i} = find(A(i, :) ~= 0);
                end
            else
                for i = 1 : num_outputs
                    indices{i} = find(A(i, 1:(i-1)) ~= 0);
                end                
            end            
        end
        
        function indices = nearestOutputIndices(this, use_outputs)
        %NEARESTOUTPUTINDICES compute nearest input (1 ... q) or output (1 ... k - 1) for each output k, k = 1 , ... , q
        % = Parameters =================================================================================================
        %   1. h (real) : stepsize
        % = Returns ====================================================================================================
        %   1. indices (vector) : array where indices(i) contains the nearest indicie. If indicies(i) < q, then the 
        %                       nearest value is the input indicies(i). If indices(i) > q, then the nearest value is the
        %                       output indicies(i) - q.
        % ==============================================================================================================
            
            num_outputs = size(this.A, 1);
            indices     = zeros(num_outputs, 1);
            z_outputs   = this.starting_times(:) + this.h;
            
            if(use_outputs)
                z_all = [this.starting_times(:); z_outputs];
                for k = 1 : num_outputs
                    indices(k) = this.nearestPoint(z_all(1 : num_outputs + k - 1), z_outputs(k));
                end
            else
                z_all = this.starting_times(:);
                for k = 1 : num_outputs
                    indices(k) = this.nearestPoint(z_all, z_outputs(k));
                end
            end

        end
        
        function index = nearestPoint(~, z_known, z_test)
        %NEARESTOUTPUTINDICES compute index of nearest value in z_known realtive to z_test
        % = Parameters =================================================================================================
        %   1. z_known (vector) - output points to consider
        %   2. z_test (bool) - point to test distance relative to
        % = Returns ====================================================================================================
        %   1. index (real) : index of nearest point
        % ==============================================================================================================
            [~, index] = min(abs(z_test - z_known));
        end

    end
    
end