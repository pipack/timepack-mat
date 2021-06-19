
% ======================================================================================================================
%
%  Abstract Block Method of the form:
%
%       y^{(n+1)} = A * y^{(n)} + r * B * y^{(n+1)} + C * y^{(n+1)} + r * D * f^{[n+1]}
%
%   where y^{(n+1)}, y^{(n)} are qx1 vectors. At the final timestep, the method produces an optional output point
%
%       y_out = a_out * y^{(n)} + r * b_out * f^{(n)} + c_out * y^{(n+1)} + r * d_out * f^{[n+1]} + r * e_out * f_out 
%
% ======================================================================================================================

classdef BLKConst < IntegratorConst
    
    properties(Abstract = true, SetAccess = protected)
        A % A Matrix
        B % B Matrix
        C % C Matrix
        D % D Matrix
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
        non_zero_A_indices = {};
        non_zero_B_indices = {};
        non_zero_C_indices = {};
        non_zero_D_indices = {};
        % -- optional output vector coefficients -----------------------------------------------------------------------
        a_out = [];
        b_out = [];
        c_out = [];
        d_out = [];
        e_out = [];
        z_out = 0;
        % -- optional coefficients for nonlinear system guess  ---------------------------------------------------------
        A_extrapolate = [];
        B_extrapolate = [];
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
        
        %Note: method is copied from IntegratorConst so that coefficient matrices can be modified
        function setClassProps(this, prop_struct) % -- allow modification of internal class properties ----------------
        
            props     = fieldnames(prop_struct);
            num_props = length(props);
            for i = 1 : num_props
                prop = props{i};
                if(ismember(prop, this.mutable_props))
                    this.(prop) = prop_struct.(prop);
                else
                    warning('property %s cannot be modified using setClassProps method.', prop);
                end                
            end

        end
        
    end
    
    methods (Access = protected) % -- generic helper functions ---------------------------------------------------------
        
        function [step_struct, y_in] = initStepStruct(this, ~, y_in, problem)
            
            num_outputs = size(this.A, 1);
            ode_dim     = size(y_in, 1);
            req_F_in    = unique(cell2mat(this.non_zero_B_indices(:).')); % all indices of f_in required for computation
            % -- all indices of f_out required for computation ---------------------------------------------------------
            req_F_out = unique(cell2mat(this.non_zero_D_indices(:).'));
            % -- all indicies of f_in required for nonlinear solve guess -----------------------------------------------
            if(this.extrapolate_initial_guess)
                req_F_extrap = find(any(this.B_extrapolate));                
                req_F_in = union(req_F_in, req_F_extrap);
            end            
            % -- flags for determining if jth RHS should be evaluated during computation -------------------------------
            req_RHS_flags = false(num_outputs, 1);
                req_RHS_flags(union(req_F_in, req_F_out)) = true;
            % -- set nearest inputs for nonlinear solve guess ----------------------------------------------------------
            nearest_output_indices = this.nearestOutputIndices( ~ this.parallel_initial_guess);
            
            step_struct = struct(                           ...
                'AT',        transpose(this.A),             ...
                'BT',        transpose(this.B),             ...
                'CT',        transpose(this.C),             ...
                'DT',        transpose(this.D),             ...
                'AT_extrp',  transpose(this.A_extrapolate), ...    
                'BT_extrp',  transpose(this.B_extrapolate), ...
                'q',         num_outputs,                   ...
                'noi',       nearest_output_indices,        ...
                'req_RHS_F', req_RHS_flags,                 ... % indices of RHS that should be evaluated during computation
                'F_in',      zeros(ode_dim, num_outputs),   ...
                'Y_out',     zeros(ode_dim, num_outputs),   ... % Stage Vector
                'F_out',     zeros(ode_dim, num_outputs)    ... % Stage Derivative Vector
            );
        
            % -- initialize F_in ---------------------------------------------------------------------------------------
            for j = (req_F_in(:)')
                conj_j = this.conjugate_inputs(j);
                if(conj_j ~= 0 && problem.real_valued) % -- check for conjugate input ----------------------------------
                    this.rhs_stats.startTimer(1); % -- add conj time on first processor
                    step_struct.F_in(:, j) = conj(step_struct.F_in(:, conj_j));
                else
                    this.rhs_stats.startTimer(j);
                    step_struct.F_in(:, j) = problem.RHS(y_in(:, j)); % -- eval directly -------------------------------
                end
                this.rhs_stats.recordRHSEval();                    
            end
            
        end     
        
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