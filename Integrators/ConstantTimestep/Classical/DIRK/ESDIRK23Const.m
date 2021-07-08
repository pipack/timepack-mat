classdef ESDIRK23Const < DI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
        gamma = (2 - sqrt(2))/2;
    end
    
    properties(SetAccess = protected)
        name  = 'ESDIRK2[3]';
        description = '2nd-order, 3-stage, L-Stable ESDIRK method (eq. 219) from C. A. Kennedy, M. H. Carpenter "Diagonally Implicit Runge-Kutta Methods for Ordinary Differential Equations. A Review," 2016'
        order = 2;
    end
    
    properties(SetAccess = protected, Dependent)
        A
        b
        c
    end
    
    methods
        
        function A = get.A(this)
            gm = this.gamma;
            b2 = (1 - 2*gm)/(4*gm);
            A  = [
                0               0       0
                gm              gm      0
                (1 - b2 - gm)   b2      gm
                ];
        end
        
        function b = get.b(this)
            gm = this.gamma;
            b2 = (1 - 2*gm)/(4*gm);
            b = [(1 - b2 - gm) b2 gm];
        end
        
        function c = get.c(this)
            c = [0 2 * this.gamma 1];
        end
        
        
        function this = ESDIRK23Const(options)
            if(nargin == 0)
                options = struct();
            end
            this@DI_RKConst(options);
        end
    end
    
end