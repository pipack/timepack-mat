% ======================================================================================================================
% Nonlinear Schrodinger Equation
%
%   iu_t + u_{xx} + u|u|^2 = 0
%
% rewritten as a real-valued system by taking u -> (u + i v) so that 
%
%   u_t = v_{xx} - (u^2 * v + v^3)
%   v_t = u_{xx} + (u^3 + v^2 * u)
%
% and with partitioning
%
%       F     : v_{xx} - (u^2 * v + v^3)
%               u_{xx} + (u^3 + v^2 * u)
%       F^{1} : v_{xx}
%               u_{xx}
%       F^{2} : -(u^2 * v + v^3)
%                (u^3 + v^2 * u)
%  
% Equation is solved in physical space using periodic boundary conditions and a spectral Fourier descritization
% ======================================================================================================================

classdef NLS_phys_rv_1d < Problem
    
	properties(SetObservable)
        tspan     = [0, 11];
        params    = struct(     ...
            'N',        2 ^ 8,  ...
            'Lx',       8 * pi ...     
        );        
    end
    
    properties(SetAccess = protected)
    	name = 'NLS';
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
    	linear_operator;
        ks; % fourier wavenumbers
    end

    methods
        
        function desc = get.description(this)
            desc = sprintf('Nonlinear Schrodinger Equation (N = %d)', this.params.N);
        end
        
        function up = RHS(this, uv, part)
        	if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % Full RHS 
                    up = this.RHS(uv, 1) + this.RHS(uv, 2);
                case 1 % Linear Term
                    Nx = this.params.N;
                    u  = uv(1 : Nx);
                    v  = uv(Nx + 1 : 2 * Nx);
                    up = [ifft((this.ks.^2) .* fft(v)); -1 * ifft((this.ks.^2) .* fft(u))];                    
                case 2 % Nonlinear Term                    
                    Nx = this.params.N;
                    u  = uv(1 : Nx);
                    v  = uv(Nx + 1 : 2 * Nx);
                    up = [-(u.^2 .* v + v.^3); (u.^3 + v.^2 .* u)]; 
            end
        end

        function J = J(this, uv, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % Full J
                    J = this.J(uv, 1) + this.J(uv, 2);
                case 1 % Linear J
                    J = this.linear_operator;
                case 2 % Nonlinear J
                    % Nonlinear Jacobian is given by:
                    %       J =  [-2 * U * V,        -3 * V^2 - U^2;
                    %            3 * U^2 + V^2,     2 * U * V       ]
                    % where U = diag(u) and V = diag(v). The following code is equivalent but faster
                    
                    Nx = this.params.N;
                    u = uv(1 : Nx);
                    v = uv(Nx + 1 : 2 * Nx);
                    
                    du = -3 * v .^ 2 - u.^2;
                    d  = -2 * u .* v;
                    dl = 3 * u .^ 2 + v.^2;
                    z  = zeros(Nx, 1);
                    
                    J = spdiags([dl d z; z -d du], [-Nx, 0, Nx], 2 * Nx, 2 * Nx);
            end
        end
        
        function Jx = Jx(this, uv, x, part)
            if(nargin == 3)
                part = 0;
            end
            
            switch part
                case 0 % Full J * x
                    Jx = this.Jx(uv, x, 1) + this.Jx(uv, x, 2);
                case 1 % Linear J * x
                    Nx = this.params.N;
                    x1 = x(1:Nx);
                    x2 = x(Nx+1:2*Nx);                    
                    Jx = [ifft((this.ks.^2) .* fft(x2)); -1 * ifft((this.ks.^2) .* fft(x1))];
                case 2
                    Nx = this.params.N;
                    x1 = x(1:Nx);
                    x2 = x(Nx+1:2*Nx);                    
                    u = uv(1 : Nx);
                    v = uv(Nx + 1 : 2 * Nx);
                    
                    du = -3 * v .^ 2 - u.^2;
                    d  = -2 * u .* v;
                    dl = 3 * u .^ 2 + v.^2;
                    
                    Jx = [d .* x1 + du .* x2; dl .*x1 - d.*x2];                    
            end
        end
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
            this.setFourierWavenumbers();
            this.setLinearOperator()
        end
        
        function setDimension(this)
            this.dimension = 2 * this.params.N;
        end
        
        function setFourierWavenumbers(this)
            N = this.params.N;                                              % num spatial points
            this.ks = [0:N/2 -N/2+1:-1]' * (2 * pi / this.params.Lx);       % Fourier wavenumbers
            if(mod(N,2) == 0)
                this.ks(N/2+1) = 0;                                         % set unbalance mode to zero
            end
        end
        
        function setInitialCondition(this)
            x  = linspace(-this.params.Lx/2, this.params.Lx/2, this.params.N + 1).'; x(end) = [];
            y0 = 1 + (1 / 100) * exp(2 * pi * 1i * x / this.params.Lx);
            this.initial_condition = [real(y0); imag(y0)];
        end

        function setLinearOperator(this)
            N  = this.params.N;
            D2 = spdiags((1i * this.ks) .^ 2, 0, N, N); % wavenumber matrix for \frac{d^2}{dx^2}
            
            DFT  = fft(eye(N));
            IDFT = ifft(eye(N));
            L    = real(IDFT * D2 * DFT);
            Z    = zeros(N);
            
            this.linear_operator = [Z -L; L Z];
        end
        
    end

end