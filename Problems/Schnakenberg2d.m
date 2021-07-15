% ======================================================================================================================
% Two-dimensional Schnakenberg Equation
%
%   u_t = gamma * (a - u + u^2 v) + (u_xx + u_yy)
%   v_t = gamma * (b - u^2 v) + d * (v_xx + v_yy)
%
% with homogeneous Neumann boundary conditions. See for example:
%
%   Eqn (26) in from Sgura, Ivonne, Benedetto Bozzini, and Deborah Lacitignola. "Numerical approximation of Turing 
%   patterns in electrodeposition by ADI methods." Journal of Computational and Applied Mathematics 236.16 (2012): 
%   4132-4147."
%
%   Luan, Vu Thai, Mayya Tokman, and Greg Rainwater. "Preconditioned implicit-exponential integrators (IMEXP) for stiff
%   PDEs." Journal of Computational Physics 335 (2017): 846-864.
%
% ======================================================================================================================

classdef Schnakenberg2d < Problem
    
    properties(SetObservable)
        tspan  = [0, 0.1];
        params = struct(        ...
            'L',        1,      ...
            'N',        64,     ...
            'a',        0.1,    ...
            'b',        0.9,    ...
            'd',        10,     ...
            'gamma',    1000    ...
        );
    end
    
    properties
        cache_reference_solutions = true;
        use_cached_reference_solutions = true;
    end
    
    properties(SetAccess = protected)
        name = 'Schnakenberg2d';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties(Access = private)
        linear_operator % stores operator for [ D_xx + D_yy , 0  ; 0 , D_xx + Dyy ]
    end
    
    methods
        
        function desc = get.description(this)
            desc = sprintf('2D Schnakenberg, N = %d^2', this.params.N);
        end
        
        function up = RHS(this, uv, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % -- full RHS -----------------------------------------------------------------------------------
                    up = this.linear_operator * uv + this.N(uv);
                case 1 % -- first component (Linear TERM) --------------------------------------------------------------
                    up = this.linear_operator * uv;
                case 2 % -- second component (Nonlinear TERM)-----------------------------------------------------------
                    up = this.N(uv);
            end
        end
        
        function J = J(this, uv, part)
            if(nargin == 2)
                part = 0;
            end
            
            if(part == 0 || part == 2)
                % -- read parameters -----------------------------------------------------------------------------------
                N = this.params.N;
                gamma = this.params.gamma;
                % -- compute nonlinear jacobian ------------------------------------------------------------------------
                N_sq = N * N;
                u = uv(1 : N_sq);
                v = uv(N_sq + 1 : end);
                u_sq = u .* u;
                u_tv = u .* v;
                Np11 = spdiags(gamma*(2 * u_tv-1), 0, N_sq,N_sq);
                Np12 = spdiags(gamma*u_sq, 0, N_sq, N_sq);
                Np21 = spdiags(gamma*(-2 * u_tv), 0, N_sq,N_sq);
                Np22 = spdiags(gamma*(-u_sq), 0, N_sq, N_sq);
                Np = [Np11, Np12; Np21, Np22];
            end
            switch part
                case 0 % -- full Jacobian ------------------------------------------------------------------------------
                    J = this.linear_operator + Np;
                case 1 % -- first component (Linear TERM) -----------------------------------------------------------
                    J = this.linear_operator;
                case 2 % -- second component (DIFFUSIVE TERM)-----------------------------------------------------------
                    J = Np;
            end
        end
        
        function N = N(this, uv)
            % -- read parameters ---------------------------------------------------------------------------------------
            Nxy = this.params.N;
            a = this.params.a;
            b = this.params.b;
            gamma = this.params.gamma;
            % -- compute nonlinearity ----------------------------------------------------------------------------------
            u = uv(1 : Nxy^2);
            v = uv(Nxy^2 + 1 : end);
            N  = [gamma*(a - u + u.^2 .* v); gamma*(b - u.^2 .* v)];
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
            x = linspace(0, this.params.L, this.params.N)';
            y = linspace(0, this.params.L, this.params.N)';
            [X, Y] = meshgrid(x, y);
            % -- steady state ------------------------------------------------------------------------------------------
            u = (this.params.a + this.params.b) * ones(this.params.N);
            v = this.params.b / (this.params.a + this.params.b)^2  * ones(this.params.N);
            % -- perturbation ------------------------------------------------------------------------------------------
            u = u + 0.0016 * cos(2*pi*(X+Y));
            v = v + 0.0016 * cos(2*pi*(X+Y));
            for k = 1:8
                u = u+0.01 * cos(2*pi*k*X);
                v = v+0.01 * cos(2*pi*k*X);
            end
            
            this.initial_condition = [u(:); v(:)];
        end
        
        function setLinearOperator(this)
            d = this.params.d;
            L = this.params.L;
            N = this.params.N;
            h = L / (N - 1);
            
            emptyMatrix = sparse(N * N, N * N);
            UxxUyy_operator = 1 / h^2 * FD_OP.UXX_UYY_N(N);
            this.linear_operator = [UxxUyy_operator, emptyMatrix; emptyMatrix, d * UxxUyy_operator];
        end
        
    end
    
end