classdef BDFConst < DI_LMMConst
    
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
        
        function this = BDFConst(options)
            if(nargin == 0)
            	options = struct();
            end
            default_field_values = {{'order', 2}};
            options = setDefaultOptions(options, default_field_values);
            this = this@DI_LMMConst(options);
            this.setOrder(options.order);       
        end
        
        function d = get.description(this)
            d = sprintf('Order %i BDF', this.order);
        end
        
        function setOrder(this, order)
            this.order            = order;
            this.name             = ['BDF', num2str(order)];
            % -- coefficients for method -------------------------------------------------------------------------------
            zs = 0 : order;
            c = double(polyCoefficients(zs, 1:order, order+1, order));
            this.a = c(1:end-1);
            this.b = [zeros(1, order) c(end)];
            % -- coefficients for nonlinear-solve guess ----------------------------------------------------------------
            zs = 0 : order - 1;
            c = double(polyCoefficients(zs, 1:order, [], order));
            this.a_extrapolate = c;
            this.b_extrapolate = zeros(1, order);
            % ----------------------------------------------------------------------------------------------------------
            this.setNonZeroVectorIndices();
            this.setStartingTimes();
        end     
        
    end

end