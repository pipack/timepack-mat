% ======================================================================================================================
% The Advection-Diffusion-Reaction (ADR) Equation
%
%   u_t = delta * u_x + epsilon * u_xx + gamma * u(u - 1/2)(1 - u)
%
% with any of the following splitings:
%
%   1. Full RHS 
%
%       F     : delta * (u_x) + epsilon * (u_xx) + gamma * u(u - 1/2)(1 - u)
%
%   2. Linear / Nonlinear (splitting = 1)
%
%           F^{1} : delta * (u_x) + epsilon * (u_xx)
%           F^{2} : gamma * u(u - 1/2)(1 - u)
%
%   3. Diffusion (splitting = 2)
%
%           F^{1} : epsilon * (u_xx)
%           F^{2} : delta * (u_x) + gamma * u(u - 1/2)(1 - u)
%
%   4. Advection (splitting = 3)
%
%           F^{1} : delta * (u_x)
%           F^{2} : epsilon * (u_xx) + gamma * u(u - 1/2)(1 - u)
%
% The equation is solved in physical space using a second-order finite differences.
%
% ======================================================================================================================

classdef ADR1d < Problem

    properties(SetObservable)
        tspan     = [0, 0.06];
        params    = struct(     ...
            'N',        100,     ...
            'epsilon',  1/100,  ...
            'delta',    -10,    ...
            'gamma',    100     ...      
        );
        splitting = 1;
    end
    
    properties
        cache_reference_solutions = true;
        use_cached_reference_solutions = true;
    end
    
	properties(SetAccess = protected)
    	name = 'ADR1d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
        xs;
    end
    
    properties(Access = private)
    	linear_operator      % stores stores epsilon * u_xx + delta * u_x
        Ux_operator          % stores delta * u_x
        Uxx_operator         % stores epsilon * u_xx
    end
        
    methods
        
        function desc = get.description(this)
            desc = sprintf('1D ADR, N = %d^2', this.params.N);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            if( part == 0 ) % -- full rhs ------------------------------------------------------------------------------ 
                up = this.RHS(u, 1) + this.RHS(u, 2);
                return;
            end
            
            switch this.splitting
                case 1 % -- semilinear with full linear operator ------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION & DIFFUSION TERM)
                            up = this.linear_operator * u;
                        case 2 % --> second component (NONLINEAR TERM)
                            up  = this.N(u);
                    end
                case 2 % -- semilinear with diffusion operator --------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (DIFFUSION TERM)
                            up = this.Uxx_operator * u;
                        case 2 % --> second component (ADVECTION NONLINEAR TERM)
                            up = this.Ux_operator * u + this.N(u);
                    end
                case 3 % -- semilinear with advection operator --------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION TERM)
                            up = this.Ux_operator * u;
                        case 2 % --> second component (DIFFUSION & NONLINEAR TERM)
                            up = this.Uxx_operator * u + this.N(u);
                    end
            end
        end

        function J = J(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            if( part == 0 ) % -- full rhs ------------------------------------------------------------------------------ 
                J = this.linear_operator + this.JN(u);
                return;
            end
            
            switch this.splitting
                case 1 % -- semilinear with full linear operator ------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION & DIFFUSION TERM)
                            J = this.linear_operator;
                        case 2 % --> second component (NONLINEAR TERM)
                            J  = this.JN(u);
                    end
                case 2 % -- semilinear with diffusion operator --------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (DIFFUSION TERM)
                            J = this.Uxx_operator;
                        case 2 % --> second component (ADVECTION NONLINEAR TERM)
                            J = this.Ux_operator + this.JN(u);
                    end
                case 3 % -- semilinear with advection operator --------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION TERM)
                            J = this.Ux_operator;
                        case 2 % --> second component (DIFFUSION & NONLINEAR TERM)
                            J = this.Uxx_operator + this.JN(u);
                    end
            end
            
        end
        
        function Jx = Jx(this, u, x, part)
            if(nargin == 3)
                part = 0;
            end
            Jx = this.J(u,part) * x;
        end
        
        function N = N(this, u)
            N = this.params.gamma * u .* (u - 1/2) .* (1 - u);
        end
        
        function JN = JN(this, u)
            d = this.params.gamma * (3*u - 1/2 - 3*u.^2);
            JN = spdiags(d, 0, this.dimension, this.dimension);
        end
        
        function L = L(this)
            L = this.linear_operator;
        end       
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
            this.setLinearOperator();
        end
        
        function setDimension(this)
            this.dimension = this.params.N;
        end
        
        function setInitialCondition(this)
            N  = this.params.N;
            x  = linspace(0, 1, N).';
            this.initial_condition = 256 * ((1 - x) .* x * 0.25).^2 + 0.3; 
            this.xs = x;
        end
        
        function setLinearOperator(this)
            epsilon = this.params.epsilon;
            delta = this.params.delta;
            N = this.params.N;
            h = 1 / (N - 1);
            
            this.Uxx_operator = epsilon / h^2 * FD_OP.UXX_N(N);
            this.Ux_operator  = delta / h * FD_OP.UX_N(N);
            this.linear_operator = this.Uxx_operator + this.Ux_operator;
        end
        
    end
    
end