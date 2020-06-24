classdef ESDIRK34Const < DI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
        gamma = 0.4358665215084589994160194;
    end
    
    properties(SetAccess = protected)
        name  = 'ESDIRK3[4]';
        description = '3rd-order 4-stage L-Stable ESDIRK method (eq. 224 from C. A. Kennedy, M. H. Carpenter "Diagonally Implicit Runge-Kutta Methods for Ordinary Differential Equations. A Review")'
        order = 3;
    end
    
    properties(SetAccess = protected, Dependent)
        A
        b
        c
    end
    
    methods
        
        function A = get.A(this)
            gm  = this.gamma;
            c3  = (3 - 20 * gm + 24 * gm^2) / (4 - 24 * gm + 24 * gm^2);
            a32 = c3 * (c3 - 2 * gm) / (4 * gm);
            b2  = (-2 + 3 * c3 + 6 * gm * (1 - c3)) / (12 * gm * (c3 - 2 * gm));
            b3  = (1 - 6 * gm + 6 * gm^2) / (3 * c3 * (c3 - 2 * gm));
            A   = [
                0                   0           0       0;
                gm                  gm          0       0;
                (c3 - a32 - gm)    a32         gm      0;
                (1 - b2 - b3 - gm)  b2          b3      gm
                ];
        end
        
        function b = get.b(this)
            gm = this.gamma;
            c3 = (3 - 20 * gm + 24 * gm^2) / (4 - 24 * gm + 24 * gm^2);
            b2 = (-2 + 3 * c3 + 6 * gm * (1 - c3)) / (12 * gm * (c3 - 2 * gm));
            b3 = (1 - 6 * gm + 6 * gm^2) / (3 * c3 * (c3 - 2 * gm));
            b  = [(1 - b2 - b3 - gm) b2 b3 gm];
        end
        
        function c = get.c(this)
            gm = this.gamma;
            c3 = (3 - 20 * gm + 24 * gm^2) / (4 - 24 * gm + 24 * gm^2);
            c  = [0 2 * this.gamma c3 1];
        end
        
        function this = ESDIRK34Const(options)
            if(nargin == 0)
                options = struct();
            end
            this@DI_RKConst(options);
        end
    end
    
end