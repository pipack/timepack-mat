classdef VanderPol < Problem
    
	properties(SetObservable)
        tspan     = [0, 0.55139];
        params    = struct(     ...
            'epsilon',  1e-1    ...    
        );        
    end
    
    properties
        cache_reference_solutions = true;
        use_cached_reference_solutions = true;
    end
        
    properties(SetAccess = protected)
    	name = 'Van der Pol';
        dimension = 2;
        initial_condition;
        real_valued = true;
        description;        
    end
    
	methods
        
        function desc = get.description(this)
            desc = sprintf('Vanderpol Oscillator (\\epsilon = %2.2g)', this.params.epsilon);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            up = [u(2); ((1 - u(1).^2).*u(2) - u(1))/this.params.epsilon];
            
            if( part == 1 )
                up(1) = 0;
            elseif ( part == 2 )
            	up(2) = 0;
            end
        end

        function J = J(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            epsilon = this.params.epsilon;
            J = [0, 1; (-2*u(1)*u(2) - 1)/epsilon, (1 - u(1)^2)/epsilon];
            
            if( part == 1 )
                J(1,2) = 0;
            elseif ( part == 2 )
            	J(2,:) = 0;
            end
        end
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
        end
        
        function setDimension(this)
            this.dimension = 2;
        end
        
        function setInitialCondition(this)
            epsilon = this.params.epsilon;
            this.initial_condition = [2, -(2 / 3) + (10 / 81) * epsilon - (292 / 2187) * epsilon^2 - (1814/19683) * epsilon^3];
        end
        
    end

end