% ======================================================================================================================
% Forced Advection-Diffusion (AD) equation
%
%   u_t = delta * (u_x + u_y) + epsilon * (u_xx + u_yy) + f(t,y)
%
% with homogenous Dirichlet boundary conditions and exact solution
%
%   u(x,y,t) = sin(2 * pi * x * (y - 1) * (1 + 2 * t)) * sin(2 * pi * y * (x - 1) * (1 + 2 * t))
%
% The following splittings are implemented:
%
%       F : delta * (u_x + u_y) + epsilon * (u_xx + u_yy) + gamma * u(u - 1/2)(1 - u)
%
%   1. Advection with forcing 
%
%       F^{1} : delta * (u_x + u_y)
%       F^{2} : epsilon * (u_xx + u_yy) + f(t,x,y)
%
%   2. Diffusion with forcing
%
%       F^{1} : epsilon * (u_xx + u_yy)
%       F^{2} : delta * (u_x + u_y) + f(t,x,y)
%
% The equation is solved in physical space using a second-order finite differences.
%
% ======================================================================================================================

classdef FADR2d < Problem

    properties(SetObservable)
        tspan  = [0, 1];
        params = struct(        ...
            'N',        32,     ...
            'epsilon',  1/100,  ...
            'delta',    -1,      ...
            'gamma',    10      ...
        );
        splitting = 1;
    end
    
	properties(SetAccess = protected)
    	name = 'FADR2d';
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
    	linear_operator         % stores stores epsilon * (u_xx + u_yy) + delta * (u_x + u_y)
        UxUy_operator           % stores delta * (u_x + u_y)
        UxxUyy_operator         % stores epsilon * (u_xx + u_yy)
        xv                      % vector of grid x-values
        yv                      % vector of grid y-values 
    end
        
    methods
        
        function desc = get.description(this)
            desc = sprintf('2D Forced AD, N = %d^2', this.params.N);
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
                case 1 % -- forcing paired with advection -------------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (DIFFUSION)
                            up = this.UxxUyy_operator * u;
                        case 2 % --> second component (ADVECTION)
                            up = this.UxUy_operator * u + this.forcingFunction(u(end)) + this.N(u);
                    end
                case 2 % -- forcing paired with diffusion -------------------------------------------------------------- 
                    switch part
                        case 1 % --> linear component (ADVECTION)
                            up = this.UxUy_operator * u;
                        case 2 % --> second component (DIFFUSION)
                            up = this.UxxUyy_operator * u + this.forcingFunction(u(end)) + this.N(u);
                    end
            end
            
        end

        function J = J(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            if( part == 0 ) % -- full rhs ------------------------------------------------------------------------------
                J = this.linear_operator + this.JForcingFunction(u(end));
                return;
            end
            
            switch this.splitting
                case 1 % -- forcing paired with advection --------------------------------------------------------------
                    switch part
                        case 1 % --> linear component (DIFFUSION)
                            J = this.UxxUyy_operator;
                        case 2 % --> second component (ADVECTION)
                            J = this.UxUy_operator + this.JforcingFunction(u(end));
                    end
                case 2 % -- forcing paired with diffusion --------------------------------------------------------------
                    switch part
                        case 1 % --> linear component (DIFFUSION)
                            J = this.UxUy_operator;
                        case 2 % --> second component (ADVECTION)
                            J = this.UxxUyy_operator + this.JforcingFunction(u(end));
                    end
            end
            
        end
        
        function Jx = Jx(this, u, x, part)
            if(nargin == 3)
                part = 0;
            end
            Jx = this.J(u,part) * x;
        end
        
        function F = forcingFunction(this, t)
            
            x = this.xv;
            y = this.yv;
            epsilon = this.params.epsilon;
            delta = this.params.delta;
            gamma = this.params.gamma;
            F = zeros(this.dimension, 1);

            F(1:end-1) = 4.*(pi+2.*pi.*t).^2.*epsilon.*(cos(2.*pi.*(1+2.*t).*(x+(-1).*y))+((-1)+ ...
            2.*x+(-2).*x.^2+2.*y+(-2).*y.^2).*cos(2.*pi.*(1+2.*t).*((-1).*y+ ...
            x.*((-1)+2.*y))))+4.*pi.*((-1)+x).*y.*cos(2.*pi.*(1+2.*t).*((-1)+ ...
            x).*y).*sin(2.*pi.*(1+2.*t).*x.*((-1)+y))+4.*pi.*x.*((-1)+y).*cos( ...
            2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y)+( ...
            -1).*gamma.*sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*(( ...
            -1)+x).*y).*(1+(-1).*sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.* ...
            pi.*(1+2.*t).*((-1)+x).*y)).*((-1/2)+sin(2.*pi.*(1+2.*t).*x.*((-1) ...
            +y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y))+(-2).*pi.*(1+2.*t).*((-1) ...
            +x+y).*delta.*sin(2.*pi.*(1+2.*t).*((-1).*y+x.*((-1)+2.*y)));
        
            F(end) = 1; % RHS for time variable
            
        end
        
        function JFF = JForcingFunction(this, t)
            
            x = this.xv;
            y = this.yv;
            epsilon = this.params.epsilon;
            delta = this.params.delta;
            gamma = this.params.gamma;
        
            data = 32.*pi.^2.*((-1)+x).*x.*((-1)+y).*y.*cos(2.*pi.*(1+2.*t).*x.*((-1) ...
            +y)).*cos(2.*pi.*(1+2.*t).*((-1)+x).*y)+(-8).*pi.^2.*(1+2.*t).*(( ...
            -1)+x+y).*((-1).*y+x.*((-1)+2.*y)).*delta.*cos(2.*pi.*(1+2.*t).*((-1) ...
            .*y+x.*((-1)+2.*y)))+16.*pi.*(pi+2.*pi.*t).*epsilon.*(cos(2.*pi.*(1+2.* ...
            t).*(x+(-1).*y))+((-1)+2.*x+(-2).*x.^2+2.*y+(-2).*y.^2).*cos(2.* ...
            pi.*(1+2.*t).*((-1).*y+x.*((-1)+2.*y))))+(-16).*pi.^2.*x.^2.*((-1) ...
            +y).^2.*sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*(( ...
            -1)+x).*y)+(-16).*pi.^2.*((-1)+x).^2.*y.^2.*sin(2.*pi.*(1+2.*t).* ...
            x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y)+(-1).*gamma.*sin(2.* ...
            pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y).*( ...
            4.*pi.*((-1)+x).*y.*cos(2.*pi.*(1+2.*t).*((-1)+x).*y).*sin(2.*pi.* ...
            (1+2.*t).*x.*((-1)+y))+4.*pi.*x.*((-1)+y).*cos(2.*pi.*(1+2.*t).* ...
            x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y)).*(1+(-1).*sin(2.* ...
            pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y))+( ...
            -1).*gamma.*sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*(( ...
            -1)+x).*y).*((-4).*pi.*((-1)+x).*y.*cos(2.*pi.*(1+2.*t).*((-1)+x) ...
            .*y).*sin(2.*pi.*(1+2.*t).*x.*((-1)+y))+(-4).*pi.*x.*((-1)+y).* ...
            cos(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).* ...
            y)).*((-1/2)+sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.* ...
            t).*((-1)+x).*y))+(-4).*pi.*((-1)+x).*y.*gamma.*cos(2.*pi.*(1+2.*t).*( ...
            (-1)+x).*y).*sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*(1+(-1).*sin(2.* ...
            pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y)).*(( ...
            -1/2)+sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*(( ...
            -1)+x).*y))+(-4).*pi.*x.*((-1)+y).*gamma.*cos(2.*pi.*(1+2.*t).*x.*(( ...
            -1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y).*(1+(-1).*sin(2.*pi.*( ...
            1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x).*y)).*((-1/2) ...
            +sin(2.*pi.*(1+2.*t).*x.*((-1)+y)).*sin(2.*pi.*(1+2.*t).*((-1)+x) ...
            .*y))+(-4).*pi.*((-1)+x+y).*delta.*sin(2.*pi.*(1+2.*t).*((-1).*y+x.*(( ...
            -1)+2.*y)))+4.*(pi+2.*pi.*t).^2.*epsilon.*((-4).*pi.*(x+(-1).*y).*sin( ...
            2.*pi.*(1+2.*t).*(x+(-1).*y))+(-4).*pi.*((-1)+2.*x+(-2).*x.^2+2.* ...
            y+(-2).*y.^2).*((-1).*y+x.*((-1)+2.*y)).*sin(2.*pi.*(1+2.*t).*(( ...
            -1).*y+x.*((-1)+2.*y))));
        
            JFF = sparse(1:this.dimension - 1, this.dimension * ones(1, this.dimension-1), data, this.dimension, this.dimension);

        end
        
        function E = exact(this, t)
        
            x = this.xv;
            y = this.yv;
            
            E = sin(2 * pi * x .* (y - 1) * (1 + 2 * t)) .* sin(2 * pi * y .* (x - 1) * (1 + 2 * t));
            
        end
        
        function N = N(this, u)
            
            N = this.params.gamma * u .* (u - 1/2) .* (1 - u);
            N(end) = 0;
            
        end
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
            this.setLinearOperator();
        end
        
        function setDimension(this)
            this.dimension = this.params.N ^ 2 + 1;
        end
        
        function setInitialCondition(this)
            
            N = this.params.N;
            xs = linspace(0, 1, N+2); xs([1 end]) = [];
            ys = linspace(0, 1, N+2); ys([1 end]) = [];
            [X, Y] = meshgrid(xs, ys);
            
            this.xv = X(:);
            this.yv = Y(:);
            
            this.initial_condition = [this.exact(0); 0];            
        end
        
        function setLinearOperator(this)
            epsilon = this.params.epsilon;
            delta = this.params.delta;
            N = this.params.N;
            h = 1 / (N + 1);
            
            this.UxxUyy_operator =  epsilon / h^2 * FD_OP.UXX_UYY_D(N);
            this.UxUy_operator   =  delta / h * FD_OP.UX_UY_D(N);
            this.linear_operator =  this.UxxUyy_operator + this.UxUy_operator;
            
            % increase size of operators by one to account for time
            t_ind = N^2 + 1;
            this.UxxUyy_operator(t_ind, t_ind) = 0;
            this.UxUy_operator(t_ind, t_ind)   = 0;
            this.linear_operator(t_ind, t_ind) = 0;
            
        end
        
    end
    
end