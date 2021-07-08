% ======================================================================================================================
% The Advection-Diffusion-Reaction (ADR) Equation
%
%   u_t = alpha * (u_x + u_y) + epsilon * (u_xx + u_yy) + gamma * u(u - 1/2)(1 - u)
%
% with homogeneous Neumann boundary conditions and any of the following splitings:
%
%   1. Full RHS 
%
%       F     : alpha * (u_x + u_y) + epsilon * (u_xx + u_yy) + gamma * u(u - 1/2)(1 - u)
%
%   2. Linear / Nonlinear (splitting = 1)
%
%           F^{1} : alpha * (u_x + u_y) + epsilon * (u_xx + u_yy)
%           F^{2} : gamma * u(u - 1/2)(1 - u)
%
%   3. Diffusion (splitting = 2)
%
%           F^{1} : epsilon * (u_xx + u_yy)
%           F^{2} : alpha * (u_x + u_y) + gamma * u(u - 1/2)(1 - u)
%
%   4. Advection (splitting = 3)
%
%           F^{1} : alpha * (u_x + u_y)
%           F^{2} : epsilon * (u_xx + u_yy) + gamma * u(u - 1/2)(1 - u)
%
% The equation is solved in physical space using a second-order finite differences.
%
% ======================================================================================================================

classdef ADR2d < Problem

    properties(SetObservable)
        tspan     = [0, 0.1];
        params    = struct(     ...
            'N',        100,     ...
            'epsilon',  1/100,  ...
            'alpha',    -10,    ...
            'gamma',    100     ...      
        );
        splitting = 1;
    end
    
	properties(SetAccess = protected)
    	name = 'ADR2d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties
        cache_reference_solutions = true;
        use_cached_reference_solutions = true;
    end
    
    properties(Access = private)
    	linear_operator         % stores stores epsilon * (u_xx + u_yy) + alpha * (u_x + u_y)
        UxUy_operator           % stores alpha * (u_x + u_y)
        UxxUyy_operator         % stores epsilon * (u_xx + u_yy)
    end
        
    methods
        
        function desc = get.description(this)
            desc = sprintf('2D ADR, N = %d^2', this.params.N);
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
                            up = this.UxxUyy_operator * u;
                        case 2 % --> second component (ADVECTION & NONLINEAR TERM)
                            up = this.UxUy_operator * u + this.N(u);
                    end
                case 3 % -- semilinear with advection operator --------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION TERM)
                            up = this.UxUy_operator * u;
                        case 2 % --> second component (DIFFUSION & NONLINEAR TERM)
                            up = this.UxxUyy_operator * u + this.N(u);
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
                            J = this.UxxUyy_operator;
                        case 2 % --> second component (ADVECTION NONLINEAR TERM)
                            J = this.UxUy_operator + this.JN(u);
                    end
                case 3 % -- semilinear with advection operator --------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION TERM)
                            J = this.UxUy_operator;
                        case 2 % --> second component (DIFFUSION & NONLINEAR TERM)
                            J = this.UxxUyy_operator + this.JN(u);
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
            this.dimension = this.params.N ^ 2;
        end
        
        function setInitialCondition(this)
            N = this.params.N;
            h = 1 / (N - 1);
            size = N * N;

            y0 = zeros(size,1);
            for j = 1:N
                y = (j-1) * h;
                for i = 1:N
                    x = (i-1) * h;
                    index = (j-1)*(N) + i;
                    y0(index) = 256 * ((1 - x) * x * (1 - y) * y)^2 + 0.3;
                end
            end
            
            this.initial_condition = y0;            
        end
        
        function setLinearOperator(this)
            epsilon = this.params.epsilon;
            alpha = this.params.alpha;
            N = this.params.N;
            h = 1 / (N - 1);
            
            this.UxxUyy_operator =  epsilon / h^2 * FD_OP.UXX_UYY_N(N);
            this.UxUy_operator   =  alpha / h * FD_OP.UX_UY_N(N);
            this.linear_operator =  this.UxxUyy_operator + this.UxUy_operator;
        end
        
    end
    
end