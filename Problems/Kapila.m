% =========================================================================
%
% Two species Kapila chemical mechanism from [1]
%
%   dy/dt = Ω                   (mole fraction Y_A)
%   dz/dt = Ω - z               (mole fraction Y_C)
%   dΘ/dt = z                   (temperature)
%
% with Ω = λyz e^{-1/(ϵΘ)} and λ = λ_0 exp(1/ϵ). Note that for numerical
% stability Ω is expressed as Ω = λ_0 yz e^{(Θ-1)/(ϵΘ)}
%
% [1] Kapila, A. K. "Homogeneous branched-chain explosion: Initiation to
%     completion." Journal of Engineering Mathematics 12.3 (1978): 221-235.
%
% =========================================================================

classdef Kapila < Problem
    
	properties(SetObservable)
        tspan  = [0, 10];
        params = struct(        ...
            'epsilon',  1e-2,   ...    
            'lambda0',   1/2     ...
        );        
    end
    
    properties
        cache_reference_solutions = true;
        use_cached_reference_solutions = true;
    end
        
    properties(SetAccess = protected)
    	name = 'Kapila';
        dimension = 3;
        initial_condition;
        real_valued = true;
        description;        
    end
    
	methods
        
        function desc = get.description(this)
            desc = sprintf('Kapila (\\epsilon = %2.2g)', this.params.epsilon);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
                       
            epsilon = this.params.epsilon;
            lambda = this.params.lambda0;            
            omega = lambda * u(1) * u(2) * exp( (u(3) - 1) / ( epsilon * u(3) ) );
            
            up = [-omega; omega - u(2); u(2)]; 
        end

        function J = J(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            epsilon = this.params.epsilon;
            lambda0 = this.params.lambda0;
            
            J = zeros(3);
            
            J(1,1) =  -lambda0 * u(2) * exp( (u(3) - 1) / ( epsilon * u(3) ) );
            J(1,2) =  -lambda0 * u(1) * exp( (u(3) - 1) / ( epsilon * u(3) ) );
            J(1,3) =  -lambda0 * u(1) * u(2) * exp( (u(3) - 1) / ( epsilon * u(3) ) ) * 1 / (epsilon * u(3)^2 );
            
            J(2,1) = - J(1,1);
            J(2,2) = - J(1,2) - 1;
            J(2,3) = - J(1,3);
            
            J(3,2) = 1;
            
        end
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
        end
        
        function setDimension(this)
            this.dimension = 3;
        end
        
        function setInitialCondition(this)
            epsilon = this.params.epsilon;
            lambda0 = this.params.lambda0;
            this.initial_condition = [1 - lambda0 * epsilon; lambda0 * epsilon; 1];
        end
        
    end

end