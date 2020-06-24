classdef ESDIRK46Const < DI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
    end
    
    properties(SetAccess = protected)
        name  = 'ESDIRK4[6]';
        description = '4th-order 6-stage stiffly-accurate L-Stable ESDIRK method (ESDIRK4(3)6L[2]SA_2 from C.A. Kennedy, M. H. Carpenter, "Diagonally implicit Runge–Kutta methods for stiff ODEs")'
        order = 4; 
        A = [
        0 0 0 0 0 0;
        31/125 31/125 0 0 0 0;
        -360286518617/7014585480527     -360286518617/7014585480527     31/125 0 0 0;
        -506388693497/5937754990171     -506388693497/5937754990171     7149918333491/13390931526268    31/125                          0                               0;
        -7628305438933/11061539393788   -7628305438933/11061539393788   21592626537567/14352247503901   11630056083252/17263101053231   31/125                          0;
        -12917657251/5222094901039      -12917657251/5222094901039      5602338284630/15643096342197    9002339615474/18125249312447    -2420307481369/24731958684496   31/125
        ];
        b = [-12917657251/5222094901039      -12917657251/5222094901039      5602338284630/15643096342197    9002339615474/18125249312447    -2420307481369/24731958684496   31/125];
        c = [0                               62/125                          486119545908/3346201505189      1043/1706                       1361/1300                       1   ];
    end
    
    methods
        
        function this = ESDIRK46Const(options)
            if(nargin == 0)
                options = struct();
            end
            this@DI_RKConst(options);
        end
    end
    
end