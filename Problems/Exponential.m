% ======================================================================================================================
% Exponential equation
%
%   u_t = \lambda y
%
% with RHS Partitioning
%
%           F     :  \lambda y
%           F^{1} :  \lambda_1 y
%           F^{2} : (\lambda - \lambda_1) y
%
% ======================================================================================================================

classdef Exponential < Problem
    
    properties(SetObservable)
        tspan     = [0, 2];
        params    = struct(     ...
            'lambda',   -1/2,   ...
            'lambda_1', 0       ...
            );
    end
    
    properties
        cache_reference_solutions = false;
        use_cached_reference_solutions = false;
    end
    
    properties(SetAccess = protected)
        name = 'Exponential';
        dimension;
        initial_condition;
        real_valued;
        description;
        
    end
    
    methods
        
        function flag = get.real_valued(this)
            if(imag(this.params.lambda) == 0)
                flag = true;
            else
                flag = false;
            end
        end
        
        function desc = get.description(this)
            desc = sprintf('Exponential, \\lambda = %f', this.params.lambda);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0
                    up = this.params.lambda * u;
                case 1
                    up = this.params.lambda_1;
                case 2
                    up = (this.params.lambda - this.params.lambda_1) * u;
                    
            end
        end
        
        function J = J(this, ~, part)
            if(nargin < 3)
                part = 0;
            end
            
            switch part
                case 0
                    J = this.params.lambda;
                case 1
                    J = this.params.lambda_1;
                case 2
                    J = (this.params.lambda - this.params.lambda_1);
                    
            end
        end
        
        function Jx = Jx(this, ~, x, part)
            if(nargin == 3)
                part = 0;
            end
            
            switch part
                case 0
                    Jx = this.params.lambda * x;
                case 1
                    Jx = this.params.lambda_1 * x;
                case 2
                    Jx = (this.params.lambda - this.params.lambda_1) * x;
            end
        end

    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
        end
        
        function setDimension(this)
            this.dimension = 1;
        end
        
        function setInitialCondition(this)
            this.initial_condition = 1;
        end
        
    end
    
end