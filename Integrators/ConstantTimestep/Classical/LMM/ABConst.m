classdef ABConst < E_LMMConst
    
    properties
        graph_line_style = {};
        eval_RHS = false;                    % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
        extrapolate_initial_guess = false;
    end
    
    properties(SetAccess = protected)
        a = [];
        b = [0];
        name  = '';
        description = '';
        order = [];
    end

    methods
        
        function this = ABConst(options)
            if(nargin == 0)
            	options = struct();
            end
            default_field_values = {{'order', 2}};
            options = setDefaultOptions(options, default_field_values);
            this = this@E_LMMConst(options);
            this.setOrder(options.order);       
        end
        
        function d = get.description(this)
            d = sprintf('Order %i Adams Bashforth', this.order);
        end
        
        function setOrder(this, order)
            this.order = order;
            this.name  = ['AB', num2str(order)];
            % -- coefficients for method -------------------------------------------------------------------------------
            zs = (-order + 1):0;
            c = double(iPolyCoefficients(zs, [0 1], 'vpa'));
            this.a = [zeros(1, order-1) 1];
            this.b = [c 0];
            % ----------------------------------------------------------------------------------------------------------
            this.setNonZeroVectorIndices();
            this.setStartingTimes();
        end     
        
    end

end