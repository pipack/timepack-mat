
% ======================================================================================================================
%
%  Diagonally Implicit - Explicit Block Method for solving
%
%       y' = fi(y) + fe(y)
%
%   The first term fe(y) will be treated implicitly while the second term fe(y) will be treated explicitly. 

%   The IMEX block method has form
%
%       y^{(n+1)} = A * y^{(n)} + r * Bi * fi^{[n]} + Be * fe^{[n]}  + C * y^{(n+1)} + r * Di * fi^{[n+1]} + r * De * fe^{[n+1]}
%
%   with optional output point 
%
%       y_out = a_out * y^{(n)} + r * bi_out * fi^{(n)} + r * be_out * fe^{(n)} + c_out * y^{(n+1)} + r * di_out * fi^{[n+1]} + r * de_out * fe^{[n+1]} + r * e_out * f_out 
%   
%  A diagonally implicit block method must satisfy the following conditions:
%   1. The matrix Di must be lower triagular
%   2. The matrix De must be strictly lower triagular
%   3. The matrix Ci must be lower triagular
%
% ======================================================================================================================

classdef DI_IMEX_BLKConst < BLKConst & IMEXConst
        
	properties(Abstract = true, SetAccess = protected)
        % method coefficient matrices
        A
        Bi
        Be
        C
        Di
        De
    end
    
    properties(Abstract = true)
        eval_RHS  % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
        linearization_index % index of input that will be used as y point for local linearization
    end
    
	properties(SetAccess = protected)
        % additional optional output coefficients
        a_out = [];
        c_out = [];
        bi_out = [];
        di_out = [];
        be_out = [];
        de_out = [];
        e_out = [];
        z_out = 0; 
        % additional properties
        non_zero_A_indices = {};
        non_zero_Bi_indices = {};
        non_zero_Be_indices = {};
        non_zero_C_indices = {};
        non_zero_Di_indices = {};
        non_zero_De_indices = {};
        % extrapolation coefficients
        A_extrapolate = [];
        Bi_extrapolate = [];
        Be_extrapolate = [];
    end
    
    methods
        
        function this = DI_IMEX_BLKConst(options)
            
            if(nargin == 0)
                options = struct();
            end
            this@BLKConst(options);
            this@IMEXConst(options);
            
            default_field_value_pairs = { ...
                {'linearization_index', 1} ...
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            this.linearization_index = options.linearization_index;
            this.verifyCefficientMatrices();
            this.setNonZeroMatrixIndices();
        end
        
        % REMARK: function body is partially copied from IntegratorConst to
        % provide access to protected variables of DI_BLKConst class. 
        % Calling the superclass method directly will not work since it 
        % cannot access protected properties (e.g. A,B,C,D matrices)
        function setClassProps(this, prop_struct)
        
            % --> copied from IntegratorConst
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
            
            % --> additional logic
            fields = fieldnames(prop_struct);
            if(contains(fields, {'A', 'Bi', 'Be', 'C', 'Di', 'De'})) 
                this.verifyCefficientMatrices();
                this.setNonZeroMatrixIndices();
            end

        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct, y_in] = initStepStruct(this, ~, y_in, problem)
            
            num_outputs = size(this.A, 1);
            ode_dim     = size(y_in, 1);
            
            [req_iRHS_flags, req_Fi_in] = this.getRequiredRHS(this.non_zero_Bi_indices, this.non_zero_Di_indices, this.Bi_extrapolate);
            [req_eRHS_flags, req_Fe_in] = this.getRequiredRHS(this.non_zero_Be_indices, this.non_zero_De_indices, this.Be_extrapolate);
                
            % -- set nearest inputs for nonlinear solve guess ----------------------------------------------------------
            nearest_output_indices = this.nearestOutputIndices( ~ this.parallel_initial_guess);
            
            step_struct = struct(                             ...
                'AT',         transpose(this.A),              ...
                'BiT',        transpose(this.Bi),             ...
                'BeT',        transpose(this.Be),             ...
                'CT',         transpose(this.C),              ...
                'DiT',        transpose(this.Di),             ...
                'DeT',        transpose(this.De),             ...
                'AT_extrp',   transpose(this.A_extrapolate),  ...    
                'BiT_extrp',  transpose(this.Bi_extrapolate), ...
                'BeT_extrp',  transpose(this.Be_extrapolate), ...
                'q',          num_outputs,                    ...
                'noi',        nearest_output_indices,         ...
                'req_iRHS_F', req_iRHS_flags,                 ... % indices of RHS that should be evaluated during computation
                'req_eRHS_F', req_eRHS_flags,                 ... % indices of RHS that should be evaluated during computation
                'Y_out',      zeros(ode_dim, num_outputs),    ... % Stage Vector
                'Fi_in',      zeros(ode_dim, num_outputs),    ... % Stage Derivative Vector
                'Fe_in',      zeros(ode_dim, num_outputs),    ... % Stage Derivative Vector
                'Fi_out',     zeros(ode_dim, num_outputs),    ... % Stage Derivative Vector
                'Fe_out',     zeros(ode_dim, num_outputs)     ... % Stage Derivative Vector
            );
        
            % -- initialize F_in ---------------------------------------------------------------------------------------
            step_struct.Fi_in(:, req_Fi_in) = this.initFin(y_in, req_Fi_in, problem, 1);
            step_struct.Fe_in(:, req_Fe_in) = this.initFin(y_in, req_Fe_in, problem, 2);
            
        end
        
        function [req_F_flags, req_F_in, req_F_out] = getRequiredRHS(this, nz_B, nz_D, B_extrap)
            num_outputs = size(this.A, 1);
            
            req_F_in  = unique(cell2mat(nz_B(:).')); % inputs
            req_F_out = unique(cell2mat(nz_D(:).')); % outputs
            if(this.extrapolate_initial_guess) 
                req_Fi_extrap = find(any(B_extrap)); % inputs for extrapolated guess               
                req_F_in = union(req_F_in, req_Fi_extrap);
            end  
            req_F = union(req_F_in, req_F_out);
            req_F_flags = false(num_outputs, 1);
            req_F_flags(req_F) = true;
        end
        
        function Fs_in = initFin(this, y_in, req_F_in, problem, part)
            
            Fs_in = zeros(size(y_in, 1), length(req_F_in));
            
            for j = (req_F_in(:)')
                conj_j = this.conjugate_inputs(j);
                if(conj_j ~= 0 && problem.real_valued) % -- check for conjugate input ----------------------------------
                    this.rhs_stats.startTimer(1); % -- add conj time on first processor
                    Fs_in(:, j) = conj(step_struct.F_in(:, conj_j));
                else
                    this.rhs_stats.startTimer(j);
                    Fs_in(:, j) = this.evalF(problem, y_in(:,j), y_in(:, this.linearization_index), part); % -- eval directly 
                end
                this.rhs_stats.recordRHSEval();                    
            end
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem)
            
            % -- read data from struct ---------------------------------------------------------------------------------          
            h     = this.h;
            q     = step_struct.q;
            AT    = step_struct.AT;
            BiT    = step_struct.BiT;
            BeT    = step_struct.BeT;
            CT    = step_struct.CT;
            DiT    = step_struct.DiT;
            DeT    = step_struct.DeT;
            AT_extrp = step_struct.AT_extrp;
            BiT_extrp = step_struct.BiT_extrp;
            BeT_extrp = step_struct.BeT_extrp;
            noi   = step_struct.noi;
            req_iRHS_flags = step_struct.req_iRHS_F;
            req_eRHS_flags = step_struct.req_eRHS_F;
            
            t_out = t_in + h;
            y_out = zeros(size(y_in));
            b_j   = zeros(size(y_in, 1), 1);
            yl    = y_in(:, this.linearization_index); % solution value for local linarization (only used if linearize = true)
            
            % -- compute outputs ---------------------------------------------------------------------------------------
            for j = 1 : q
            	
                this.step_stats.startTimer(j);
                
                if(problem.real_valued == false || this.conjugate_outputs(j) == 0)
                    this.linear_solver.stats.setIndex(j);
                        
                    % -- b_j = y_in * AT(:, j) + h * fe_in * BeT(:, j) + h * fi_in * BiT(:, j)+ y_out * CT(:, j) + fe_out * DeT(:, 1:j-1) + fi_out * DiT(:, 1:j-1) ---------
                    b_j(:) = 0;    
                    if(~isempty(this.non_zero_A_indices{j}))
                        b_j = b_j + y_in(:, this.non_zero_A_indices{j}) * AT(this.non_zero_A_indices{j}, j);
                    end
                    if(~isempty(this.non_zero_Bi_indices{j}))
                        b_j = b_j + h * (step_struct.Fi_in(:, this.non_zero_Bi_indices{j}) * BiT(this.non_zero_Bi_indices{j}, j));
                    end
                    if(~isempty(this.non_zero_Be_indices{j}))
                        b_j = b_j + h * (step_struct.Fe_in(:, this.non_zero_Be_indices{j}) * BeT(this.non_zero_Be_indices{j}, j));
                    end
                    if(~isempty(this.non_zero_C_indices{j}))
                        b_j = b_j + y_out(:, this.non_zero_C_indices{j}) * CT(this.non_zero_C_indices{j}, j);
                    end
                    if(~isempty(this.non_zero_Di_indices{j}))
                        b_j = b_j + h * (step_struct.Fi_out(:, this.non_zero_Di_indices{j}))  * DiT(this.non_zero_Di_indices{j}, j);
                    end
                    if(~isempty(this.non_zero_De_indices{j}))
                        b_j = b_j + h * (step_struct.Fe_out(:, this.non_zero_De_indices{j}))  * DeT(this.non_zero_De_indices{j}, j);
                    end
                    
                    c_j = h * this.Di(j, j);
                    
                    if(c_j ~= 0) % -- solve linear system y_j = b_j + c_j * L(y_j) ----------------------------------
                        
                        % -- choose guess for nonlinear system ---------------------------------------------------------
                        if(this.extrapolate_initial_guess)
                            y_guess = zeros(size(y_in, 1), 1);
                            if(~isempty(this.A_extrapolate))
                                y_guess = y_guess + y_in * AT_extrp(:, j);
                            end
                            if(~isempty(this.Bi_extrapolate))
                                y_guess = y_guess + h * step_struct.Fi_in * BiT_extrp(:, j);
                            end
                            if(~isempty(this.Be_extrapolate))
                                y_guess = y_guess + h * step_struct.Fe_in * BeT_extrp(:, j);
                            end
                        else % choose nearest value as guess
                            if(noi(j) <= q)
                                y_guess = y_in(:, noi(j));
                            else
                                y_guess = y_out(:, noi(j) - q);
                            end
                        end
                        
                        % -- solve system: Note: first ODE component is assumbed to be linear, so problem.J returns
                        % linear operator regardlesss of expansion point y_in(:,1)
                        if(problem.real_valued && this.real_valued_outputs(j)) % -- output is real valued time point -------------
                        	[y_out(:,j), clean_exit] = this.solveBC(problem,real(b_j), real(c_j), y_guess, yl);
                        else % -- output is complex valued -------------------------------------------------------------
                        	[y_out(:,j), clean_exit] = this.solveBC(problem,b_j, c_j, y_guess, yl);
                        end

                        emergency_exit = ((~clean_exit) || any(isinf(y_out(:,j))) || any(isnan(y_out(:,j))));
                        if(emergency_exit)
                            y_out = NaN;
                            return;
                        end

                    else % -- system is explicit -----------------------------------------------------------------------
                        if(problem.real_valued && this.real_valued_outputs(j))
                            y_out(:,j) = real(b_j);
                        else
                            y_out(:,j) = b_j;
                        end
                    end
                
                else
                    conj_j = this.conjugate_outputs(j);
                    y_out(:, j) = conj(y_out(:, conj_j));
                end
                    
                % -- eval RHS ------------------------------------------------------------------------------------------
                if(req_iRHS_flags(j))
                    step_struct.Fi_out(:, j) = this.evalFiOut(j, y_out(:,j), step_struct.Fi_out, problem, yl, b_j, c_j);
                end
                if(req_eRHS_flags(j))
                    step_struct.Fe_out(:, j) = this.evalFeOut(j, y_out(:,j), step_struct.Fe_out, problem, yl);
                end
                
                this.rhs_stats.recordRHSEval();
                this.step_stats.recordStep();
            end
            
            % -- update F_in for next step -----------------------------------------------------------------------------
            if(this.linearize) % re-evaluate RHS using jacobian evaluated at output (instead of input)
                yl = y_out(:, this.linearization_index);
                for j = 1 : q
                    if(req_iRHS_flags(j))
                        step_struct.Fi_in(:, j) = this.evalFiOut(j, y_out(:,j), step_struct.Fi_out, problem, yl);
                    end
                    if(req_eRHS_flags(j))
                        step_struct.Fe_in(:, j) = this.evalFeOut(j, y_out(:,j), step_struct.Fe_out, problem, yl);
                    end                    
                end
            else
                step_struct.Fi_in = step_struct.Fi_out;
                step_struct.Fe_in = step_struct.Fe_out;
            end
            
        end
         
        function F_out = evalFeOut(this, j, y_out, y_out_prev, problem, yl)
            
            this.rhs_stats.startTimer(j);
            conj_j = this.conjugate_outputs(j);
            if(conj_j ~= 0 && problem.real_valued) % -- check for conjugate input ------------------------------
                F_out = conj(y_out_prev(:, conj_j));
            else
                F_out = this.evalF(problem, y_out, yl, 2); %problem.RHS(y_out, 2);
            end
        end
        
        function F_out = evalFiOut(this, j, y_out, y_out_prev, problem, yl, b_j, c_j)
            
            this.rhs_stats.startTimer(j);
            conj_j = this.conjugate_outputs(j);
            if(conj_j ~= 0 && problem.real_valued) % -- check for conjugate input ------------------------------
                F_out = conj(y_out_prev(:, conj_j));
            else
                if(this.eval_RHS || nargin == 6 || c_j == 0)
                    F_out = this.evalF(problem, y_out, yl, 1); % -- eval directly ------- % problem.RHS(y_out, 1); 
                else % -- compute via re-arrangement -----------------------------------------------------------
                    if(problem.real_valued && this.real_valued_outputs(j))
                        F_out = real(y_out - b_j) / real(c_j); 
                    else
                        F_out = (y_out - b_j) / c_j; 
                    end
                end
            end
            
        end
        
        function [t_user, y_user] = userOutput(this, t_in, y_in, step_struct_in, t_out, y_out, step_struct_out, problem)
            
            flag = ~isempty(this.a_out) && ~isempty(this.bi_out) && ~isempty(this.be_out) && ~isempty(this.c_out) && ~isempty(this.di_out) && ~isempty(this.de_out) && ~isempty(this.e_out);
            if(flag)
                
                % -- read data from struct -----------------------------------------------------------------------------
                h = this.h;
                q = step_struct_out.q;
                
                req_iRHS_flags = step_struct_out.req_iRHS_F;
                req_eRHS_flags = step_struct_out.req_eRHS_F;
                yl = y_in(:, this.linearization_index);
                
                % -- check if F input requires additional Evaluations (Implicit component) -----------------------------
                uncomputed_Fi_in  = setdiff(find(this.bi_out), find(req_iRHS_flags));
                for j = uncomputed_Fi_in(:)'
                    this.rhs_stats.startTimer(j);
                    step_struct_out.Fi_in(:, j) = this.evalF(problem, y_in(:,j), yl, 1); % problem.RHS(y_in(:, j), 1);
                    this.rhs_stats.recordRHSEval();
                end
                % -- check if F output requires additional Evaluations (Implicit component) ----------------------------
                uncomputed_Fi_out = setdiff(find(this.di_out), find(req_iRHS_flags));
                for j = uncomputed_Fi_out(:)'
                    this.rhs_stats.startTimer(j);
                    step_struct_out.Fi_out(:, j) = this.evalF(problem, y_out(:,j), yl, 1); %problem.RHS(y_out(:, j), 1);
                    this.rhs_stats.recordRHSEval();
                end
                % -- check if F input requires additional Evaluations --------------------------------------------------
                uncomputed_Fe_in  = setdiff(find(this.be_out), find(req_eRHS_flags));
                for j = uncomputed_Fe_in(:)'
                    this.rhs_stats.startTimer(j);
                    step_struct_out.Fi_in(:, j) = this.evalF(problem, y_in(:,j), yl, 2); %problem.RHS(y_in(:, j), 2);
                    this.rhs_stats.recordRHSEval();
                end
                % -- check if F output requires additional Evaluations -------------------------------------------------
                uncomputed_Fe_out = setdiff(find(this.de_out), find(req_eRHS_flags));
                for j = uncomputed_Fe_out(:)'
                    this.rhs_stats.startTimer(j);
                    step_struct_out.Fi_out(:, j) = this.evalF(problem, y_out(:,j), yl, 2); %problem.RHS(y_out(:, j), 2);
                    this.rhs_stats.recordRHSEval();
                end
                
                % -- form nonlinear equation y = b + c * f(y) --------------------------------------------------------------
                b = y_in * this.a_out(:) + h * (step_struct_out.Fi_in * this.bi_out(:) + step_struct_out.Fe_in * this.be_out(:)) + ...
                    y_out * this.c_out(:) + h * (step_struct_out.Fi_out * this.di_out(:) + step_struct_out.Fe_out * this.de_out(:));
                c = h * this.e_out;
                
                if(c ~= 0) % -- solve linear system ---------------------------------------------------------------------
                    
                    % -- choose guess for nonlinear system -----------------------------------------------------------------
                    nearest_index = this.nearestPoint([this.starting_times(:); this.starting_times(:) + this.h], this.z_out + this.h);
                    if(nearest_index <= q)
                        y_guess = y_in(:, nearest_index);
                    else
                        y_guess = y_out(:, nearest_index - q);
                    end
                    % -- solve system --------------------------------------------------------------------------------------
                    [y_user, clean_exit] = this.linear_solver.solveBC(problem.J(y_in(:,1),1), b, c, y_guess);
                    if(imag(this.z_out) == 0 && problem.real_valued)
                        y_user = real(y_user);
                    end
                    emergency_exit = ((~clean_exit) || any(isinf(y_user)) || any(isnan(y_user)));
                    if(emergency_exit)
                        t_user = NaN;
                        y_user = NaN;
                        return;
                    end
                else % -- system is explicit -------------------------------------------------------------------------------
                    if(imag(this.z_out) == 0 && problem.real_valued)
                        y_user = real(b);
                    else
                        y_user = b;
                    end
                end
                
                t_user = problem.tspan(end) + h * this.z_out; % output value
                
            else
                t_user = t_in + this.h;
                y_user = y_in;
            end
            
        end
        
    end
    
    methods (Access = protected)
                
        function verifyCefficientMatrices(this)
            isstril = @(A) istril(A) && ~any(diag(A)); % tests if A is strictly lower triangular
            if(~istril(this.De))
                error('Diagonally Implicit IMEX Block Method must have a lower triangular Di Matrix');
            end
            if(~isstril(this.Di))
                error('Diagonally Implicit IMEX Block Method must have a strictly lower triangular De Matrix');
            end
            if(~isstril(this.C))
                error('Diagonally Implicit IMEX Block Method must have a strictly lower triangular C Matrix');
            end
        end
        
        function setNonZeroMatrixIndices(this)
        %SETNONZEROSTAGEINDICES set nonzero row entries of the block matrices A, B, C, D
        % = Parameters =================================================================================================
        % = Returns ====================================================================================================
        % ==============================================================================================================   
            this.non_zero_A_indices  = this.nonzeroMatrixRows(this.A, true);
            this.non_zero_Bi_indices = this.nonzeroMatrixRows(this.Bi, true);
            this.non_zero_Be_indices = this.nonzeroMatrixRows(this.Be, true);
            this.non_zero_C_indices  = this.nonzeroMatrixRows(this.C, false);
            this.non_zero_Di_indices = this.nonzeroMatrixRows(this.Di, false);
            this.non_zero_De_indices = this.nonzeroMatrixRows(this.De, false);
        end

    end
    
end