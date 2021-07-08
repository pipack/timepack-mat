classdef IMEXRK324L2SA < DI_IMEXRKConst
    
    properties
        graph_line_style = {};
        eval_RHS = true;
        linearly_implicit = false;
        linearize = false;
    end
    
    properties(SetAccess = protected)
        name  = 'ARK3(2)4L[2]SA';
        description = '3rd-order 4-stage IMEX-RK method from C. A. Kennedy and M. H. Carpenter, Additive Runge-Kutta schemes for convection-diffusion-reaction equations, Appl. Numer. Math., 44 (2003), pp. 139-181.'
        order = 3;
    end
    
    properties(SetAccess = protected)
        Ai = [
            0,                                          0,                                       0,                                           0;                                   
            1767732205903.0 / 4055673282236.0,          1767732205903.0 / 4055673282236.0,       0,                                           0;
            2746238789719.0 / 10658868560708.0,      - 640167445237.0 / 6845629431997.0,         1767732205903.0 / 4055673282236.0,           0;
            1471266399579.0 / 7840856788654.0,       - 4482444167858.0 / 7529755066697.0,        11266239266428.0 / 11593286722821.0,         1767732205903.0 / 4055673282236.0;
        ]
        Ae = [
            0,                                          0,                                        0,                                        0;
            1767732205903.0 / 2027836641118.0,          0,                                        0,                                        0;
            5535828885825.0 / 10492691773637.0,         788022342437.0 / 10882634858940.0,        0,                                        0;
            6485989280629.0  / 16251701735622.        - 4246266847089.0 / 9704473918619.0,        10755448449292.0 / 10357097424841,        0;
        ]
        bi = [1471266399579.0 / 7840856788654.0, - 4482444167858.0 / 7529755066697.0, 11266239266428.0 / 11593286722821.0, 1767732205903.0 / 4055673282236.0]
        be = [1471266399579.0 / 7840856788654.0, - 4482444167858.0 / 7529755066697.0, 11266239266428.0 / 11593286722821.0, 1767732205903.0 / 4055673282236.0]
        c  = [ 0, 1767732205903 / 2027836641118, 3 / 5, 1 ];
    end
    
    methods
        
        function this = IMEXRK324L2SA(options)
            if(nargin == 0)
                options = struct();
            end
            options = setDefaultOptions(options, {{'linearize', false}});
            this@DI_IMEXRKConst(options);
            this.linearize = options.linearize;
        end
    end

end