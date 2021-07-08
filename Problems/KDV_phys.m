% ======================================================================================================================
% KDV Equation
%
%   u_t + \delta^2 * u_{xxx} + \frac{1}{2}(u^2)_x = 0
%
% with partitioning
%
%       F     : -\delta^2 * u_{xxx} - \frac{1}{2}(u^2)_x
%       F^{1} : -\delta^2 * u_{xxx}
%       F^{2} : - \frac{1}{2}(u^2)_x
%
% The equation is solved in physical space using periodic boundary conditions and a spectral Fourier discretization.
% ======================================================================================================================

classdef KDV_phys < Problem
    
    properties(SetObservable)
        tspan     = [0, 3.6/pi];
        params    = struct(     ...
            'N',        2 ^ 8,  ...
            'delta',    0.022,  ...
            'Lx',       2       ...
            );
    end
    
    properties
        cache_reference_solutions = true;
        use_cached_reference_solutions = true;
    end
    
    properties(SetAccess = protected)
        name = 'KDV';
        dimension;
        initial_condition;
        real_valued = true;
        description;
    end
    
    properties(Access = private)
        linear_operator;        % full linear operator
        ks;                     % fourier wave numbers
        DDX                     % linear operator for D/DX
    end
    
    methods
        
        function desc = get.description(this)
            desc = sprintf('KDV (N = %d)', this.params.N);
        end
        
        function up = RHS(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % Full RHS
                    up = this.RHS(u, 1) + this.RHS(u, 2);
                case 1 % Linear Term
                    up = -1 * (this.params.delta ^ 2) * ifft( -1i  * ( this.ks ) .^ 3 .* fft(u));
                case 2 % Nonlinear Term
                    up = (-1/2) * ifft( 1i * this.ks .* fft((u.^ 2)));
            end           
        end
        
        function J = J(this, u, part)
            if(nargin == 2)
                part = 0;
            end
            
            switch part
                case 0 % Full J
                    J = this.J(u, 1) + this.J(u, 2);
                case 1 % Linear Term J
                    J = this.linear_operator;
                case 2 % Nonlinear Term J
                    N = this.params.N;
                    J = - this.DDX * spdiags(u, 0, N, N);
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
                    Jx = -1 * (this.params.delta ^ 2) * real( ifft( -1i  * ( this.ks ) .^ 3 .* fft(x)));
                case 2 % Nonlinear Term Jx
                    Jx = -1 * real(ifft(1i * this.ks .* fft(u.* x)));
            end
        end
        
    end
    
    methods(Access = protected)
        
        function reset(this, varargin)
            this.setInitialCondition();
            this.setFourierWavenumbers();
            this.setLinearOperator();
        end
        
        function setDimension(this)
            this.dimension = this.params.N;
        end
        
        function setFourierWavenumbers(this)
            N = this.params.N;                                              % num spatial points
            this.ks = [0:N/2 -N/2+1:-1]' * (2 * pi / this.params.Lx);       % Fourier wavenumbers
            if(mod(N,2) == 0)
                this.ks(N/2+1) = 0;                                         % set unbalance mode to zero
            end        
        end
        
        function setInitialCondition(this)
            x  = linspace(0, this.params.Lx, this.params.N + 1).'; x(end) = [];
            this.initial_condition = cos(pi * x);
        end
        
        function setLinearOperator(this) % used for forming full Jacobian matrix
            
            N  = this.params.N;
            D  = spdiags((1i * this.ks), 0, N, N);                                % wavenumber matrix for \frac{d}{dx}
            D3 = spdiags((1i * this.ks) .^ 3, 0, N, N);                           % wavenumber matrix for \frac{d^3}{dx^3}
            
            DFT  = fft(eye(N));
            IDFT = ifft(eye(N));
            
            this.linear_operator = -this.params.delta^2 * real(IDFT * D3 * DFT);
            this.DDX             = real(IDFT * D * DFT);
            
        end
                
    end
    
end