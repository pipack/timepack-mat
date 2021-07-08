classdef IMEXRK1 < DI_IMEXRKConst
    
    properties
        graph_line_style = {};
        eval_RHS = true;
        linearly_implicit = false;
        linearize = false;
    end
    
    properties(SetAccess = protected)
        name  = 'IMEX-RK1';
        description = 'Forward-backward Euler'
        order = 1;
        Ai = [0 0; 0 1]
        bi = [0 1];
        Ae = [0 0; 1 0];
        be = [1 0];
        c  = [0 1];
    end
    
    methods
          
        function this = IMEXRK1(options)
            if(nargin == 0)
                options = struct();
            end
            options = setDefaultOptions(options, {{'linearize', false}});
            this@DI_IMEXRKConst(options);
            this.linearize = options.linearize;
        end
    end
    
end