% ======================================================================================================================
% Burgers' Equation
%
%   u_t = \nu * u_xx + (1/2) * (u^2)_x
%
% With RHS Partitioning
%
%           F     : \nu * u_xx + (1/2) * (u^2)_x
%           F^{1} : \nu * u_xx
%           F^{2} : (1/2) * (u^2)_x
%
%   Equation is solved in physical space using 2nd order finite difference discretization
%
% ======================================================================================================================

classdef Burgers1d < Problem
    
    properties(SetObservable)
        tspan     = [0, 1];
        params    = struct(     ...
            'N',        2000,   ...
            'Lx',       1,      ...
            'nu',       3e-4    ...
            );
    end
    
    properties(SetAccess = protected)
        name = 'Burgers1d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties(Access = private)
        linear_operator;
    end
    
    methods
        
        function desc = get.description(this)
            desc = sprintf('Viscous Burgers 1d, \nu = %2.2g, N = %d', this.params.nu, this.params.N);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % -- full RHS -----------------------------------------------------------------------------------
                    up = this.RHS(u, 1) + this.RHS(u, 2);
                case 1 % -- linear component (ADVECTION & DIFFUSION TERM) ----------------------------------------------
                    up = this.linear_operator * u;
                case 2 % -- second component (NONLINEAR TERM) ----------------------------------------------------------
                    N    = this.params.N;
                    h    = this.params.Lx / (N + 1);
                    c    = -1 / (4 * h);
                    temp = ones(N, 1);
                    M    = c * (spdiags(temp, 1, N, N) - spdiags(temp, -1, N, N));
                    up   = M * (u.^2);
            end            
        end
        
        function J = J(this, u, part)
        	if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % -- full RHS -----------------------------------------------------------------------------------
                    J = this.J(u, 1) + this.J(u, 2);
                case 1 % -- linear component (ADVECTION & DIFFUSION TERM) ----------------------------------------------
                    J = this.linear_operator;
                case 2 % -- second component (NONLINEAR TERM) ----------------------------------------------------------
                    N    = this.params.N;
                    h    = this.params.Lx / (N + 1);
                    c    = -1 / (2 * h);
                    J    = c * (spdiags(u, 1, N, N) - spdiags(u, -1, N, N));
            end 
        end
        
        function Jx = Jx(this, u, x, part)
            if(nargin == 3)
                part = 0;
            end
            
            switch part
                case 0 % Full Jx
                    Jx = this.Jx(u, x, 1) + this.Jx(u, x, 2);
                case 1 % Linear Term Jx
                    Jx = this.linear_operator * x;
                case 2 % Nonlinear Term Jx
                    N    = this.params.N;
                    h    = this.params.Lx / (N + 1);
                    c    = -1 / (2 * h);
                    Jx = c * (spdiags(u, 1, N, N) * x - spdiags(u, -1, N, N) * x);
            end
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
            N = this.params.N;
            x = linspace(0, this.params.Lx, this.params.N + 2);
            y0 = zeros(N, 1);
            for i = 2 : N + 1
                y0(i-1,1) = (sin(3 * pi * x(i))) .^ 2 .* (1 - x(i)) .^ (3/2);
            end
            this.initial_condition = y0;
        end
        
        function setLinearOperator(this)
            N = this.params.N;
            h = this.params.Lx / (N + 1);
            c = this.params.nu / h ^ 2;
            this.linear_operator = c * FD_OP.UXX_D(N);
        end
        
    end
    
end