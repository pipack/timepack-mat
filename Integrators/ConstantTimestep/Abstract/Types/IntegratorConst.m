% ======================================================================================================================
% INTEGRATORCONST : The base class for all integrators with constant stepsize. It automatically computes initial 
%                   conditions for methods and performs the integration by calling the step function.
%
% -- Public Properties -------------------------------------------------------------------------------------------------
%
% 1. num_timesteps  : numer of timesteps for integrating problem
%
% 2. starting_times : vector of starting times relative in natural coordinates tau = (t - t_0)/h.
%                    
% 3. starting_integrator : starting integrator for computing initial times
%
% 4. starting_integrator_timeout : maximum time for determining initial conditions. 
%                                 
% 5. max_num_store : integer cooresponding to the maximum number of solution points to store; default is one.
%
% 6. ind_store : specific indices to store (overrides num_store if specified ).
%
% -- Public Methods ----------------------------------------------------------------------------------------------------
%
%  1. solve(problem) - integrates problem object
%  2. resetStats     - resents all integrator stats object 
%  3. getStats       - returns struct containing all stats objects 
%
% ======================================================================================================================

classdef IntegratorConst < handle
    
    properties (Abstract = true)
        graph_line_style
    end
    
    properties (Abstract = true, SetAccess = protected)
        name
        description
        starting_times      % array of required starting times
        order               % method order
    end
    
    properties (Access = protected)
        h
        mutable_props % determine which class properties can be changed using setClassProps method
    end
    
    properties
        % -- integrator properties -------------------------------------------------------------------------------------
        num_timesteps = 1;
        starting_integrator;
        starting_integrator_timeout; % max time in seconds for computing each initial condition
        % -- solution storage properties -------------------------------------------------------------------------------
        max_num_store = 1;
        ind_store = [];
        % -- stats objects ---------------------------------------------------------------------------------------------
        step_stats
        rhs_stats
    end
    
    methods (Abstract = true, Access = protected)
        initStepStruct(this, t_in, y_in, problem);                              % initializes data for timestepping
        step(this, t_in, y_in, step_struct, problem);                           % advances timestep
        userOutput(t_in, y_in, struct_in, t_out, y_out, struct_out, problem);   % returns user output
    end
    
    methods
        
        function this = IntegratorConst(options)
            
            default_field_value_pairs = { ...
                {'starting_integrator', []} ...
                {'starting_integrator_timeout', 600} ...
                {'graph_line_style', {'k.--'}} ...
                {'record_stats',  false} ...
                {'max_num_store', 1} ...
                {'ind_store', []} ...
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            
            % -- set props ---------------------------------------------------------------------------------------------
            props = {'starting_integrator', 'starting_integrator_timeout', 'graph_line_style'};
            this.setEmptyClassProps(props, options);
            
            % -- initialize stats object -------------------------------------------------------------------------------
            this.step_stats = TimestepStats();
                this.step_stats.record = options.record_stats;            
            this.rhs_stats = RHSStats();
                this.rhs_stats.record = options.record_stats;
        end
        
        function set.num_timesteps(this, num_timesteps)
            %SET.NUM_TIMESTEP ensure that num_timesteps is valid
            input_is_int = floor(num_timesteps) == num_timesteps;
            intput_is_positive = num_timesteps > 0;
            if(input_is_int && intput_is_positive)
                this.num_timesteps = num_timesteps;
            else
                warning('Integrator: num_timesteps must be strictly positive integer');
            end
        end
        
        function set.ind_store(this, ind_store)
            %SET.NUM_TIMESTEP ensure that num_timesteps is valid
            input_is_int = all(floor(ind_store) == ind_store);
            input_is_positive = all(ind_store > 1);
            if(input_is_int && input_is_positive)
                this.ind_store = ind_store;
            else
                warning('Integrator: ind_num_store must empty of array of integer smaller than num_steps');
            end
        end
        
        function [t_user, y_user] = solve(this, problem)
            if(isempty(this.num_timesteps))
                error('Integrator: number of timesteps has not been set.');
            end
            % -- initialize output variables ---------------------------------------------------------------------------
            [y_user, t_user, i_store, store_index] = initOutputVars(this, problem);            
            % -- set stepsize ------------------------------------------------------------------------------------------
            this.setStepsize(problem);
            % -- obtain initial conditions -----------------------------------------------------------------------------
            [t_out, y_out, clean_exit] = this.initialConditions(problem);
            if(~clean_exit)
                [t_user(store_index:end), y_user(:, store_index:end)] = deal(NaN);
                warning('Integrator: Emergency exit during initial condition computation.');
                return;
            end
            % -- solve problem -----------------------------------------------------------------------------------------
            struct_out = this.initStepStruct(t_out, y_out, problem);
            num_steps = this.num_timesteps;
            for i = 1 : num_steps
                % --> step
                [t_in, struct_in, y_in] = deal(t_out, struct_out, y_out);
                [t_out, y_out, struct_out] = this.step(t_in, y_in, struct_in, problem);
                % --> emergency exit conditions
                if(any(isinf(y_out(:))) || any(isnan(y_out(:))))
                    [t_user(store_index:end), y_user(:, store_index:end)] = deal(NaN);
                    warning('Integrator: Emergency exit at step %i', i);
                    break;
                end
                % --> store output
                if(i == i_store(store_index))
                    [t_user(store_index), y_user(:, store_index)] = this.userOutput(t_in, y_in, struct_in, t_out, y_out, struct_out, problem);
                    store_index = store_index + 1;
                end                
            end
            % -- output solution if only one output argument -----------------------------------------------------------
            if(nargout <= 1)
                t_user = y_user;
            end
        end
        
        function resetStats(this)
            this.step_stats.reset();
            this.rhs_stats.reset();
            if(isprop(this, 'nonlinear_solver'))
                this.nonlinear_solver.stats.reset();
                this.nonlinear_solver.linear_solver.stats.reset();
            end
            if(isprop(this, 'linear_solver'))
                this.linear_solver.stats.reset();
            end
            if(isprop(this, 'phi_evaluator'))
                this.phi_evaluator.stats.reset();
            end
        end
        
        function stats_struct = getStats(this)
            stats_struct = struct( ...
                'timestep', copy(this.step_stats), ...
                'rhs', copy(this.rhs_stats) ...
                );
            if(isprop(this, 'nonlinear_solver'))
            	stats_struct.nonlinear = {copy(this.nonlinear_solver.stats), copy(this.nonlinear_solver.linear_solver.stats)};
            end
            if(isprop(this, 'linear_solver'))
                stats_struct.linear = copy(this.linear_solver.stats);
            end
            if(isprop(this, 'phi_evaluator'))
                stats_struct.phi = copy(this.phi_evaluator.stats);
            end
        end
        
        function setClassProps(this, prop_struct) % -- allow modification of internal method properties ----------------
        
            props     = fieldnames(prop_struct);
            num_props = length(props);
            for i = 1 : num_props
                prop = props{i};
                if(ismember(prop, this.mutable_props))
                    this.(prop) = prop_struct.(prop);
                else
                    warning('property %s cannot be modified using setClassProps method.');
                end                
            end

        end
        
    end
    
    methods(Access = protected)
        
        function setStepsize(this, problem)
            this.h = (problem.tspan(end) - problem.tspan(1))/this.num_timesteps;
        end
        
        % == START FUNCTIONS FOR COMPUTING INITIAL INPUTS ==============================================================
        
        function [t_initial, y_initial, clean_exit] = initialConditions(this, problem, scaling_factor)
            %INITIALCONDITIONS Summary of this function goes here
            %   Detailed explanation goes here
            if(nargin == 2)
                scaling_factor = this.h;
            end
            
            % define local timeout function for ode15s
            function status = timeout(start_time, max_duration)
                if(toc(start_time) > max_duration)
                   error('timeout:Interupt', 'Computation exceeded allowed time interval!') 
                end
                status = 0; % function must return 0 for integrator to progress
            end
            
            clean_exit  = true;
            t_initial   = problem.tspan(1) + scaling_factor * [0; this.starting_times(:)]; % prepent initial condition \tau = 0
            num_outputs = length(t_initial) - 1;
            ode_dim     = length(problem.initial_condition);
            y_initial   = zeros(ode_dim, num_outputs+1); % prenent initial condition
            y_initial(:,1) = problem.initial_condition;
            
            integration_order = startingValueOrder(this, this.starting_times);
            
            for i=1:num_outputs
                % -- start and end integration times -----------
                start_index = integration_order(1,i) + 1;
                end_index   = integration_order(2,i) + 1;
                t_start     = t_initial(start_index);
                y_start     = y_initial(:, start_index);
                t_end       = t_initial(end_index);
                % -- wrap into mutable problem ---------------------
                problem_i = this.initCondProblemWrapper(problem, t_start, t_end, y_start);
                if(diff(problem_i.tspan) == 0)
                    y_end = y_start;
                else
                    if(~isempty(this.starting_integrator))
                        % -- set starting integrator num_timesteps ---
                        if(~isempty(this.order))
                            est_error = max(1e-16, (this.h ^ this.order)/100); % assume local error of integrator at each step is proprotional to h^order / 100
                        else
                            est_error = 1e-16;
                        end                        
                        this.starting_integrator.num_timesteps = max(10, ceil(diff(problem_i.tspan)/(est_error ^ (1/this.starting_integrator.order)))); % increased min to 10 steps for ADR2D
                        y_end = this.starting_integrator.solve(problem_i);
                    else % USE ODE15
                        % -- set up timeout function -------------------------------------------------------------------
                        start_time = tic();
                        outputFun = @(t, y, flag) timeout(start_time, this.starting_integrator_timeout);
                        % -- compute initial condition -----------------------------------------------------------------
                        try
                            ode_options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14, 'Jacobian', @(varargin) problem_i.J(varargin{2:end}), 'OutputFcn', outputFun);
                            [~, sol] = ode15s(@(varargin) problem_i.RHS(varargin{2:end}), problem_i.tspan, problem_i.initial_condition, ode_options);
                            y_end = sol(end,:);
                        catch E
                            if strcmp(E.identifier,'timeout:Interupt')
                                warning('Initial condition computation timed out');
                                y_end = NaN;
                            else
                            	rethrow(E);
                            end
                        end
                    end
                end
                % -- emergency exit condition --------------------------------------------------------------------------
                if(any(isinf(y_end(:))) || any(isnan(y_end(:))))
                    clean_exit = false;
                    break;
                end
                y_initial(:, end_index) = y_end;
            end
            % remove \tau = 0 initial condition
            y_initial = y_initial(:,2:end);
            t_initial = t_initial(2:end);
        end
        
        function integration_indices = startingValueOrder(~, times)
            % OPTIMALSTARTINGINTEGRATIONINDICES  computes optimal order to compute output_times. Used by initialConditions
            % PARAMETERS
            %   times (vector) : vector containing output times in local coordinates.
            % RETURNS
            % integration_indices (array) : returns array of size 2 x length(output_times) indicating optimal way to
            %                               compute values times: times(integration_indices(i,2)) should be computed
            %                               using times(integration_indices(i,1)) as initial condition value for i = 1, ...,
            %                               length(output_times). If integration_indices(i,1) = 0, then t = 0 should be
            %                               initial condition.
            
            num_outputs = length(times);
            
            % -- ic solve order --
            z_known = zeros(1, num_outputs + 1); % known z values
            z_known_index = zeros(1, num_outputs + 1); % index of known z values
            
            z_unknown = double(times); %
            z_unknown_index = 1:num_outputs;
            integration_indices = zeros(2, num_outputs);
            
            for i=1:num_outputs
                % -- find closest z_known to any z_unknown ---------
                min_dist = Inf;
                min_start_zindex = -1;
                min_end_zindex = -1;
                min_end_index = -1;
                for j=1:i
                    [min_dist_cand, index_cand] = min(abs(z_known(j) - z_unknown)); % closest z_unknown to z_known(j)
                    if(min_dist_cand < min_dist)
                        min_dist         = min_dist_cand;
                        min_start_zindex = z_known_index(j);
                        min_end_zindex   = z_unknown_index(index_cand);
                        min_end_index    = index_cand;
                    end
                end
                % -- move closest z_unknown to z_known ---
                z_known(i + 1) = times(min_end_zindex);
                z_known_index(i + 1) = min_end_zindex;
                % -- remove unknown point ----------------
                z_unknown(min_end_index) = [];
                z_unknown_index(min_end_index) = [];
                % -- store start and end index -----------
                integration_indices(:,i) = [min_start_zindex, min_end_zindex];
            end
            
        end
        
        function [problem_ic] = initCondProblemWrapper(~, problem, t_start, t_end, y_start)
            %COMPLEXPROBLEMWRAPPER converts a problem with complex time integration bounds into a initial value problem with real
            % integation bounds. This is accomplished by wrapping the original problem into a MutableProblem object.
            
            % -- compute real integration path -------------------------------------------------------------------------
            z           = t_end - t_start;
            r           = abs(z);
            exp_i_theta = z/abs(z);
            tspan       = [0 r];
            % -- wrap problem ------------------------------------------------------------------------------------------
            problem_ic = MutableProblem(struct('tspan', tspan, 'dimension', length(y_start), 'initial_condition', y_start, 'real_valued', problem.real_valued));
            if(problem.has('RHS'))
                problem_ic.RHS = @(varargin) exp_i_theta * problem.RHS(varargin{:});
            end
            if(problem.has('J'))
                problem_ic.J = @(varargin) exp_i_theta * problem.J(varargin{:});
            end
            if(problem.has('Jx'))
                problem_ic.Jx =  @(varargin) exp_i_theta * problem.Jx(varargin{:});
            end
        end
        
        function [y_user, t_user, i_store, store_index] = initOutputVars(this, problem)
            
            % --> get solution indices that should be stored
            i_store = this.getIndStore(this.num_timesteps, this.max_num_store, this.ind_store);
            store_index = 1;
            % --> initialize storage arrays            
            y_user = zeros(problem.dimension, length(i_store));
            t_user = zeros(length(i_store), 1);
            % --> if needed, store initial condition           
            if(i_store(1) == 0)
                t_user(store_index)   = problem.tspan(1);
                y_user(:,store_index) = problem.initial_condition;
                store_index = 2;
            end            
        
        end
               
        function inds = getIndStore(~, num_steps, max_store, ind_store)
            % SETINDSTORE determines the indices of the timesteps that 
            % should should be stored and returned to the user
            
            max_store = min(num_steps + 1, max_store);
            
            if(~isempty(ind_store)) % case 1 : user has specified ind_store
                inds = sort(ind_store, 'asc');
                if(inds(end) ~= num_steps)
                    warning('final solution will not be returned. Check value of integrator ind_store parameter');
                end
            elseif(max_store == 1) % case 2 : num_store = 1 (only store last timestep)    
                inds = num_steps;
            else % case 2 : user has specified num_store > 1
                delta = floor(num_steps / (max_store - 1)); % divisor
                shift = mod(num_steps, max_store - 1);      % remainder
                skips = [0, delta * ones(1, max_store - 1)] + [0, ones(1, shift), zeros(1, max_store - 1 - shift)]; % distribute remainder over initial steps
                inds  = cumsum(skips);
            end
        end
        
        function setEmptyClassProps(this, props, values)
            for i = 1 : length(props)
                if(isempty(this.(props{i})))
                    this.(props{i}) = values.(props{i});
                end
            end
        end
        
    end
    
end