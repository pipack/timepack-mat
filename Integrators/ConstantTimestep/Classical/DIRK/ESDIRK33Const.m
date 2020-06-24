classdef ESDIRK33Const < DI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
        gamma = (3 + sqrt(3))/6;
        c3    = 1;
    end
    
    properties(SetAccess = protected)
        name  = 'ESDIRK3[3]';
        description = '3rd-order 3-stage A-Stable ESDIRK method (eq. 222 from C. A. Kennedy, M. H. Carpenter "Diagonally Implicit Runge-Kutta Methods for Ordinary Differential Equations. A Review")'
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
            c3_ = this.c3;
            a32 = c3_ * (c3_ - 2 * gm) / (4 * gm);
            
            A  = [
                0                   0       0
                gm                  gm      0
                (c3_ - a32 - gm)     a32     gm
                ];
        end
        
        function b = get.b(this)
            gm = this.gamma;
            c3_ = this.c3;
            b2  = (-2 + 3 * c3_) / (12 * (c3_ - 2 * gm) * gm);
            b3  = (1 - 3 * gm) / (3 * c3_ * (c3_ - 2 * gm));
            b = [(1 - b2 - b3) b2 b3];
        end
        
        function c = get.c(this)
            c3_  = this.c3;
            c = [0 2 * this.gamma c3_];
        end
        
        
        function this = ESDIRK33Const(options)
            if(nargin == 0)
                options = struct();
            end
            this@DI_RKConst(options);
        end
    end
    
end