classdef IMEXRK436L2SA < DI_IMEXRKConst
    
    properties
        graph_line_style = {};
        eval_RHS = true;
        linearly_implicit = false;
        linearize = false;
    end
    
    properties(SetAccess = protected)
        name  = 'ARK4(3)6L[2]SA';
        description = '4rd-order 6-stage IMEX-RK method from C. A. Kennedy and M. H. Carpenter, Additive Runge-Kutta schemes for convection-diffusion-reaction equations, Appl. Numer. Math., 44 (2003), pp. 139-181.'
        order = 4;
    end
    
    properties(SetAccess = protected)
        Ai = [
            0,                                    0,                              0,                               0,                            0,                      0;
            0.25,                                 0.25,                           0,                               0,                            0,                      0;
            8611.0 / 62500.0,                   - 1743.0 / 31250.0,               0.25,                            0,                            0,                      0;
            5012029.0 / 34652500.0,             - 654441.0 / 2922500.0,           174375.0 / 388108.0,             0.25,                         0,                      0;
            15267082809.0 / 155376265600.0,     - 71443401.0 / 120774400.0,       730878875.0 / 902184768.0,       2285395.0 / 8070912.0,        0.25,                   0;
            82889.0 / 524892.0                    0,                              15625.0 / 83664.0,               69875.0 / 102672.0,         - 2260.0 / 8211.0,        0.25              
        ];
        Ae = [
            0,                                    0,                                        0,                                         0,                                         0,                       0;       
            0.5,                                  0,                                        0,                                         0,                                         0,                       0;
            13861.0 / 62500.0,                    6889.0 / 62500.0,                         0,                                         0,                                         0,                       0;
          - 116923316275.0 / 2393684061468.0,   - 2731218467317.0 / 15368042101831.0,       9408046702089.0 / 11113171139209.0,        0,                                         0,                       0;
          - 451086348788.0 / 2902428689909.0,   - 2682348792572.0 / 7519795681897.0,        12662868775082.0 / 11960479115383.0,       3355817975965.0 / 11060851509271.0,        0,                       0;
            647845179188.0 / 3216320057751.0,     73281519250.0 / 8382639484533.0,          552539513391.0 / 3454668386233.0,          3354512671639.0 / 8306763924573.0,         4040.0 / 17871.0         0;
        ]
        bi = [82889.0 / 524892.0, 0.0, 15625.0 /  83664.0, 69875.0 / 102672.0, - 2260.0 / 8211.0, 0.25];
        be = [82889.0 / 524892.0, 0.0, 15625.0 /  83664.0, 69875.0 / 102672.0, - 2260.0 / 8211.0, 0.25];
        c  = [ 0, 0.5, 83 / 250, 31 / 50, 17 / 20, 1 ];
    end
    
    methods
        
        function this = IMEXRK436L2SA(options)
            if(nargin == 0)
                options = struct();
            end
            options = setDefaultOptions(options, {{'linearize', false}});
            this@DI_IMEXRKConst(options);
            this.linearize = options.linearize;
        end
    end
    
end