classdef Problem < matlab.mixin.Copyable % (See https://www.mathworks.com/help/matlab/ref/matlab.mixin.copyable-class.html)
    
	properties (Abstract = true, SetObservable)
        tspan;
        params;        
    end
    
    properties (Abstract = true)
        cache_reference_solutions
        use_cached_reference_solutions
    end
    
    properties(SetAccess = protected)
    	reference_solution = [];
    end
    
	properties(Abstract, SetAccess = protected)
    	dimension;
        name;
        description;
        initial_condition;
        real_valued;
    end

    methods(Abstract = true, Access = protected)
        reset(this);
        setDimension(this);
    end
    
    methods
        
        function this = Problem()
            addlistener(this,'tspan','PostSet',@this.ObservablesHaveChanged);
            addlistener(this,'params','PostSet',@this.ObservablesHaveChanged);
            this.setDimension();
            this.reset();
        end
        
        function flag = has(this, name)
            if(ismethod(this, name))
                flag = true;
            elseif(isprop(this, name) && ~isempty(this.(name)))
                flag = true;
            else
                flag = false;
            end
        end
        
        function sol = get.reference_solution(this)
            if(isempty(this.reference_solution))
                this.setRefSolution();
            end            
            sol = this.reference_solution;
        end
        
        function Jx = Jx(this, u, x, part) % default Finite difference implementation in case no jacobian is provided
            if(nargin == 3)
                part = 0;
            end
            epsilon = sqrt(eps) / norm(x, 2); % Q.N. Ludlow, D.K. Shaw, 'A matrix-free preconditioned Newton/GMRES method for unsteady Navier-Stokes solutions'
            Jx = (this.RHS(u + epsilon * x, part) - this.RHS(u, part)) / epsilon;
        end
        
    end
    
    methods(Access = protected)
    
        function ObservablesHaveChanged(this, varargin)
            this.reference_solution = [];
            this.setDimension();
            this.reset();
        end
        
        function setRefSolution(this)
            cache_path = this.refSolutionCachePath();
            % -- load from cache ---------------------------------------------------------------------------------------
            if(this.use_cached_reference_solutions && isfile(cache_path))
                load(cache_path, 'reference_solution');
                this.reference_solution = reference_solution;
                return
            end
            
            if(isfield(this.params, 'ref_solve_tol'))
                ref_solve_tol = this.params.ref_solve_tol;
            else
                ref_solve_tol = 3e-14;
            end
            options = odeset('Jacobian', @(varargin) this.J(varargin{2:end}), 'RelTol', ref_solve_tol, 'AbsTol', ref_solve_tol);
            timeSpan = [this.tspan(1), (this.tspan(1) + this.tspan(2))/2, this.tspan(2)];
            [~,exsol] = ode15s(@(varargin) this.RHS(varargin{2:end}), timeSpan, this.initial_condition, options);
            rs = transpose(exsol(end,:));           
            this.reference_solution = rs;
            
            % -- cache reference solution ------------------------------------------------------------------------------
            if(this.cache_reference_solutions)
                [~, ~] = mkdir(fileparts(cache_path)); % ensure directory exists;
                reference_solution = rs;
                save(cache_path, 'reference_solution');
            end
            
        end
        
        function hash = identifier(this)
            s = this.params;
            s.tspan = this.tspan;
            hash = MD5DataHash(s);
        end
        
        function cache_path = refSolutionCachePath(this)
            pipack_problems_dir = fileparts(which('Problem'));
            cache_path = fullfile(pipack_problems_dir, 'cache', this.name, [this.identifier(), '.mat']);
        end

    end
    
end