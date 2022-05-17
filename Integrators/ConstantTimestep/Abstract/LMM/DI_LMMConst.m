
% ======================================================================================================================
%  Abstract Class for Linear Multistep Methods: 
%
%       y_{n+1} = \sum_{j=1}^p a(j) y_{n-p+j} + h * \sum_{j=1}^{p+1} b(j) * f_{n-p+j}
%
%   Properties
%       a (vector) - solution coefficients:   a(j) * y_{n-p+j}   for j = 1 ... p
%       b (vector) - derivative coefficients: h b(j) * f_{n-p+j} for j = 1 ... p+1
%
%       Examples:
%           Implicit Euler  ->  a = [1],     b = [0 1]
%           AB2             ->  a = [1],     b = [1/2 1/2]
%           AB3             ->  a = [0 1],   b = [-1/12 2/3 5/12]
% ======================================================================================================================

classdef DI_LMMConst < IntegratorConst & ImplicitIntegratorConst
    
    properties(Abstract = true, SetAccess = protected)
        a % vector of coeffients for y
        b % vector of coefficients for f
    end
    
    properties(Abstract = true)
        eval_RHS                    % boolean - if true RHS will evaluated directly, otherwise F will obtained algebraically after nonlinear solve
        extrapolate_initial_guess   % if true coefficients a_extrapolate and b_extrapolate will be used to form initial guess for any nonlinear systems 
    end
    
    properties(SetAccess = protected)
        starting_times = 0;
        non_zero_a_indices = [];
        non_zero_b_indices = [];
        % -- optional coefficients for nonlinear system guess  ---------------------------------------------------------
        a_extrapolate = [];
        b_extrapolate = [];
    end
    
    properties(Access = protected)
        step_backwards_for_ic = false;   % if true, initial conditions will be computed behind first timestep
    end
    
    methods
        
        function this = DI_LMMConst(options) % -- Constructor -------------------------------------------------------------
            this = this@IntegratorConst(options);
            this = this@ImplicitIntegratorConst(options);
            this.verifyCoefficientMatrices();
            this.setNonZeroVectorIndices();
            this.setStartingTimes();
        end
        
        function set.step_backwards_for_ic(this, val)
            this.step_backwards_for_ic = val;
            this.setStartingTimes();
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct] = initStepStruct(this, ~, y_in, problem)
            
            p           = length(this.starting_times);
            ode_dim     = size(y_in, 1);
            b_prev      = this.b(1:end-1);
            compute_rhs = ~isempty(this.non_zero_b_indices) || (this.extrapolate_initial_guess && sum(this.b_extrapolate ~= 0) > 0);
            
            step_struct = struct(                                    ...
                'a',            this.a(this.non_zero_a_indices),     ...
                'b_out',        this.b(end),                         ... % coefficient for f_{n+1}
                'b_prev',       b_prev(this.non_zero_b_indices),     ... % coefficients for f_{n-j} for j = 0, ..., p.
                'F',            zeros(ode_dim, p),                   ...
                'compute_rhs',  compute_rhs,                         ...
                'ode_dim',      ode_dim                              ...
            );
        
            step_struct.a = step_struct.a(:); % ensure column vector
            step_struct.b_prev = step_struct.b_prev(:); % ensure column vector
            
            % -- initialize F_in ---------------------------------------------------------------------------------------
            if(compute_rhs)
                for j = 1 : p
                    this.rhs_stats.startTimer()
                    step_struct.F(:, j) = problem.RHS(y_in(:, j));
                    this.rhs_stats.recordRHSEval();
                end
            end
            
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem)
            
            % -- read data from struct ---------------------------------------------------------------------------------
            h         = this.h;
            ac        = step_struct.a;
            bc_prev   = step_struct.b_prev;
            bc_out    = step_struct.b_out;
            ode_dim   = step_struct.ode_dim;
            
            t_out = t_in + h;
            y_out = zeros(size(y_in));
            
            this.step_stats.startTimer();
            
            % -- compute X ---------------------------------------------------------------------------------------------
            X = zeros(ode_dim, 1);
            if(~isempty(this.non_zero_a_indices))
                X = X + y_in(:, this.non_zero_a_indices) * ac;
            end
            if(~isempty(this.non_zero_b_indices))
                X = X + step_struct.F(:, this.non_zero_b_indices) * (h * bc_prev);
            end
            
            % -- solve nonlinear system --------------------------------------------------------------------------------
            if(bc_out ~= 0) % implicit                
                if(this.extrapolate_initial_guess) % choose guess for nonlinear system                 
                    ae = this.a_extrapolate(:);
                    be = this.b_extrapolate(:);
                    y_guess = y_in * ae + step_struct.F * (h * be);
                else % choose nearest value as guess
                    y_guess = y_in(:, end);
                end
                [y_out(:, end), clean_exit] = this.nonlinear_solver.solveBC(problem, X, h * bc_out, y_guess);
            else % explicit
                y_out(:, end) = X;
                clean_exit = true;
            end
            y_out(:, 1:end-1) = y_in(:, 2:end);
            
            
            emergency_exit = ((~clean_exit) || any(isinf(y_out(:,end))) || any(isnan(y_out(:, end))));
            if(emergency_exit)
                y_out = NaN;
                return;
            end
            
            % -- eval RHS ----------------------------------------------------------------------------------------------
            if(step_struct.compute_rhs)
                step_struct.F(:, 1:end-1) = step_struct.F(:, 2:end);
                if(this.eval_RHS || bc_out == 0)
                    this.rhs_stats.startTimer()
                    step_struct.F(:, end) = problem.RHS(y_out(:, end)); % -- eval directly ---------------------------------
                    this.rhs_stats.recordRHSEval();
                else
                    step_struct.F(:, end) = (y_out(:,end) - X) / (h * bc_out); % -- compute via re-arrangement -------------
                end
            end
            
            this.step_stats.recordStep();         
        end
        
        function [t_user, y_user] = userOutput(this, t_in, y_in, struct_in, t_out, y_out, struct_out, problem)
            t_user = t_out(end);
            y_user = y_out(:, end);
        end
        
    end
    
    methods (Access = protected)
        
        function setStepsize(this, problem)
            p = max(length(this.a), length(this.b) - 1);
            if(this.step_backwards_for_ic)
                this.h = (problem.tspan(end) - problem.tspan(1))/(this.num_timesteps);
            else
                this.h = (problem.tspan(end) - problem.tspan(1))/(p - 1 + this.num_timesteps);
            end
        end
        
        function setNonZeroVectorIndices(this)
            tol = eps * 10;
            this.non_zero_a_indices = find(abs(this.a) > tol);
            this.non_zero_b_indices = find(abs(this.b(1:end-1)) > tol); % non-zero coefficients for f_{n-j} j = 0, 1, ... p
        end
        
        function setStartingTimes(this)
            p = max(length(this.a), length(this.b) - 1);
            if(this.step_backwards_for_ic)
                this.starting_times = -p+1:0;
            else
                this.starting_times = 0:p-1;
            end            
        end
        
        function verifyCoefficientMatrices(this)
            if(length(this.a) ~= length(this.b) - 1)
                error('Length of coefficient vector "b" (cooresponding to derivatives) must be one larger than length of coefficient vector "a" (coorespondng to values)');
            end
        end
                
    end
    
end
