classdef AMConst < DI_LMMConst
    
    properties
        graph_line_style = {};
        eval_RHS = false;                    % boolean - if true RHS will evaluated directly, otherwise F will algebraically obtained after nonlinear solve
        extrapolate_initial_guess = false;
    end
    
    properties(SetAccess = protected)
        a = [1]
        b = [0 1]
        name  = '';
        description = '';
        order = [];
    end

    methods
        
        function this = AMConst(options)
            if(nargin == 0)
            	options = struct();
            end
            default_field_values = {{'order', 2}};
            options = setDefaultOptions(options, default_field_values);
            this = this@DI_LMMConst(options);
            this.setOrder(options.order);       
        end
        
        function d = get.description(this)
            d = sprintf('Order %i Adams Moulton', this.order);
        end
        
        function setOrder(this, order)
            this.order = order;
            this.name  = ['AM', num2str(order)];
            % -- coefficients for method -------------------------------------------------------------------------------
            zs = (-order + 2):1;
            c = double(iPolyCoefficients(zs, [0 1], 'vpa'));
            this.a = [zeros(1, order-2) 1];
            this.b = c;
            % -- coefficients for nonlinear-solve guess ----------------------------------------------------------------
            zs = (-order + 2):0;
            c = double(polyCoefficients(zs, 1:(order-1), [], 1));
            this.a_extrapolate = c;
            this.b_extrapolate = zeros(1, order - 1);
            % ----------------------------------------------------------------------------------------------------------
            this.setNonZeroVectorIndices();
            this.setStartingTimes();
        end     
        
    end

end