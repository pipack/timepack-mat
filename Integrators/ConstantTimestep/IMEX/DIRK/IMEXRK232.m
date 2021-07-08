
classdef IMEXRK232 < DI_IMEXRKConst
    
    properties
        graph_line_style = {};
        eval_RHS = true;
        gamma = (2 - sqrt(2))/2;
        delta = -2*sqrt(2)/3;
        linearly_implicit = false;
        linearize = false;
    end
    
    properties(SetAccess = protected)
        name  = 'IMEX-RK(2,3,2)';
        description = '2nd-order 3-stage IMEX-RK method from U. M. Ascher, S. J. Ruuth, and R. J. Spiteri, Implicit-Explicit Runge-Kutta methods for time-dependent partial differential equations, Appl. Numer. Math., 25 (1997), pp. 151-167.'
        order = 2;
    end
    
    properties(SetAccess = protected, Dependent)
        Ai
        Ae
        bi
        be
        c
    end
    
    methods
        
        function Ai = get.Ai(this)
            gm = this.gamma;
            Ai = [
                0    0        0;
                0    gm       0;
                0    1 - gm   gm
            ];
        end
        
        function bi = get.bi(this)
            gm = this.gamma;
            bi = [0, 1 - gm, gm];
        end
        
        function Ae = get.Ae(this)
            gm = this.gamma;
            dl = this.delta;
            Ae = [
                0     0        0;
                gm    0        0;
                dl    1 - dl   0
            ];
        end
        
        function be = get.be(this)
            gm   = this.gamma;
            be = [0, 1 - gm, gm];
        end
        
        function c = get.c(this)
            c    = [0 this.gamma 1];
        end
        
        function this = IMEXRK232(options)
            if(nargin == 0)
                options = struct();
            end
            options = setDefaultOptions(options, {{'linearize', false}});
            this@DI_IMEXRKConst(options);
            this.linearize = options.linearize;
        end
    end

end