classdef ESDIRK45Const < DI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
        gamma = 0.5728160624821348554080014;
        c3 = 1308256777188 / 2690004194437;
        c4 = 2026389075477 / 2726940318254;
    end
    
    properties(SetAccess = protected)
        name  = 'ESDIRK4[5]';
        description = '4th-order 5-stage L-Stable ESDIRK method (eq. 233 and eq. (241)-(246) from C. A. Kennedy, M. H. Carpenter "Diagonally Implicit Runge-Kutta Methods for Ordinary Differential Equations. A Review")'
        order = 4;      
    end
    
    properties(SetAccess = protected, Dependent)
        A
        b
        c
    end
    
    properties(Access = private)
        b2  = @(gm, c3, c4) (3 - 12 * gm + 4 * c4 * (-1 + 3*gm) - 2 * c3 * (2 - 6 * gm + c4 * (-3 + 6 * gm))) / (24 * gm * (2 * gm - c3) * (2 * gm - c4));
        b3  = @(gm, c3, c4, ps1, ps2) (ps2 - 4 * c4 * ps1) / (12 * c3 * (c3 - c4) * (c3 - 2 * gm));
        b4  = @(gm, c3, c4, ps1, ps2) (ps2 - 4 * c3 * ps1) / (12 * c4 * (c4 - c3) * (c4 - 2 * gm));
        a32 = @(gm, c3) (c3 * (c3 - 2 * gm)) / (4 * gm);
        a42 = @(gm, c3, c4, ps1, ps2, ps3, ps4) (c4 * (c4 - 2 * gm) * (-4 * c3^2 * ps1 - 2 * gm * ps2 + c3 * ps3 + 2 * c4 * ps4 )) / (4 * gm * (2 * gm - c3) * (4 * c3 * ps1 - ps2));
        a43 = @(gm, c3, c4, ps1, ps2, ps4) ((c4 - c3) * c4 * (c4 - 2 * gm) * ps4) / (c3 * (c3 - 2 * gm) * (4 * c3 * ps1 - ps2));
        ps1 = @(gm) (1 - 6 * gm + 6 * gm^2)
        ps2 = @(gm) (3 - 20 * gm + 24 * gm^2)
        ps3 = @(gm) (5 - 36 * gm + 48 * gm^2)
        ps4 = @(gm) (-1 + 12 * gm - 36 * gm^2 + 24 * gm^3);        
    end
    
    methods(Access = private)
        function [b2, b3, b4, a32, a42, a43, ps1, ps2, ps3, ps4] = coeff(this, gamma_, c3_, c4_)
        	ps1 = this.ps1(gamma_);
            ps2 = this.ps2(gamma_);
            ps3 = this.ps3(gamma_);
            ps4 = this.ps4(gamma_);
            b2  = this.b2(gamma_, c3_, c4_);
            b3  = this.b3(gamma_, c3_, c4_, ps1, ps2);
            b4  = this.b4(gamma_, c3_, c4_, ps1, ps2);
            a32 = this.a32(gamma_, c3_);
            a42 = this.a42(gamma_, c3_, c4_, ps1, ps2, ps3, ps4);
            a43 = this.a43(gamma_, c3_, c4_, ps1, ps2, ps4);
        end
        
    end
    
    methods
        
        function A = get.A(this)
            gm  = this.gamma;
            c3_ = this.c3;
            c4_ = this.c4;
            [b2, b3, b4, a32, a42, a43] = this.coeff(gm, c3_, c4_);

            A   = [
                0                       0       0       0      0;   
                gm                      gm      0       0      0;
                (c3_ - a32 - gm)        a32     gm      0      0;
                (c4_ - a42 - a43 - gm)  a42     a43     gm     0;
                (1 - b2 - b3 - b4 - gm) b2      b3      b4     gm;
                ];
        end
        
        function b = get.b(this)
            gm  = this.gamma;
            c3_ = this.c3;
            c4_ = this.c4;
            [b2, b3, b4] = this.coeff(gm, c3_, c4_);            
            
            b  = [(1 - b2 - b3 - b4 - gm), b2, b3, b4, gm];
        end
        
        function c = get.c(this)
            gm  = this.gamma;
            c3_ = this.c3;
            c4_ = this.c4;
            
            c  = [0 2*gm c3_ c4_ 1];
        end
        
        function this = ESDIRK45Const(options)
            if(nargin == 0)
                options = struct();
            end
            this@DI_RKConst(options);
        end
    end
    
end