classdef IMEXRK427L2SA < DI_IMEXRKConst
    
    properties
        graph_line_style = {};
        eval_RHS = true;
        linearly_implicit = false;
        linearize = false;
    end
    
    properties(SetAccess = protected)
        name  = 'ARK4(3)7L[2]SA';
        description = '4rd-order 7-stage IMEX-RK method from Kennedy, Christopher A., and Mark H. Carpenter. "Higher-order additive Runge–Kutta schemes for ordinary differential equations." Applied Numerical Mathematics 136 (2019): 183-205.'
        order = 4;
    end
    
    properties(SetAccess = protected)
        Ai = [
            0,                                     0,                                     0,                                     0,                                     0,                                     0,                                    0;
            1235/10000,                            1235/10000,                            0,                                     0,                                     0,                                     0,                                    0;
            624185399699 / 4186980696204,          624185399699 / 4186980696204,          1235/10000,                            0,                                     0,                                     0,                                    0;
            1258591069120 / 10082082980243,        1258591069120 / 10082082980243,      - 322722984531 / 8455138723562,          1235/10000,                            0,                                     0,                                    0;
          - 436103496990 / 5971407786587,        - 436103496990 / 5971407786587,        - 2689175662187 / 11046760208243,        4431412449334 / 12995360898505,        1235/10000,                            0,                                    0;
          - 2207373168298 / 14430576638973,      - 2207373168298 / 14430576638973,        242511121179 / 3358618340039,          3145666661981 / 7780404714551,         5882073923981 / 14490790706663,        1235/10000,                           0;
            0,                                     0,                                     9164257142617 / 17756377923965,      - 10812980402763 / 74029279521829,       1335994250573 / 5691609445217,         2273837961795 / 8368240463276,        1235/10000;  
        ];
        Ae = [
            0,                                     0,                                     0,                                     0,                                    0,                                    0,                                     0;
            247 / 1000,                            0,                                     0,                                     0,                                    0,                                    0,                                     0;
            247 / 4000,                            2694949928731 / 7487940209513,         0,                                     0,                                    0,                                    0,                                     0;
            464650059369 / 8764239774964,          878889893998 / 2444806327765,        - 952945855348  / 12294611323341,        0,                                    0,                                    0,                                     0;
            476636172619 / 8159180917465,        - 1271469283451 / 7793814740893,       - 859560642026  / 4356155882851,         1723805262919 / 4571918432560,        0,                                    0,                                     0;
            6338158500785 / 11769362343261,      - 4970555480458 / 10924838743837,        3326578051521 / 2647936831840,       - 880713585975 / 1841400956686,       - 1428733748635 / 8843423958496,        0,                                     0;
            760814592956 / 3276306540349,          760814592956 / 3276306540349,        - 47223648122716 / 6934462133451,        71187472546993 / 9669769126921,     - 13330509492149 / 9695768672337,       11565764226357 / 8513123442827,        0;
        ]
        bi = [0, 0, 9164257142617 / 17756377923965, - 10812980402763 / 74029279521829, 1335994250573 / 5691609445217, 2273837961795 / 8368240463276, 1235/10000 ];
        be = [0, 0, 9164257142617 / 17756377923965, - 10812980402763 / 74029279521829, 1335994250573 / 5691609445217, 2273837961795 / 8368240463276, 1235/10000 ];
        c  = [ 0, 247 / 1000, 4276536705230 / 10142255878289, 67 / 200, 3 / 40, 7 / 10, 1 ];
    end
    
    methods
        
        function this = IMEXRK427L2SA(options)
            if(nargin == 0)
                options = struct();
            end
            options = setDefaultOptions(options, {{'linearize', false}});
            this@DI_IMEXRKConst(options);
            this.linearize = options.linearize;
        end
    end
       
end