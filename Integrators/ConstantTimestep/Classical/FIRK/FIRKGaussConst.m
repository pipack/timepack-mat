classdef FIRKGaussConst < FI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
        extrapolate_initial_guess = false;
    end
    
    properties(SetObservable)
        s = 2; % number of stages (user can dynamically set)
    end
    
    properties(SetAccess = protected, Dependent)
        A
        b
        c
        A_extrapolate
        name
        description
        order
    end
        
    methods
        
        function this = FIRKGaussConst(options)
            if(nargin == 0)
                options = struct();
            end
            default_options = {{'s', 2}};
            options = setDefaultOptions(options, default_options);            
            this@FI_RKConst(options);
            
            addlistener(this, 's', 'PostSet', @this.ObservablesHaveChanged);            
            this.s = options.s;
        end
        
        function A = get.A(this)
            z = legpts(this.s, [0,1]);
            ab = [zeros(this.s, 1), z(:)];
            A = double(iPolyCoefficients(z, ab, 'vpa'));
        end
        
        function b = get.b(this)
            z = legpts(this.s, [0,1]);
            ab = [0, 1];
            b = double(iPolyCoefficients(z, ab, 'vpa'));
        end
        
        function c = get.c(this)
            c = legpts(this.s, [0,1]);
        end
        
        function Ae = get.A_extrapolate(this)
            z = legpts(this.s, [0,1]);
            z_aug = [0; z; 1]; % input, stage, output
            Ae = double(polyCoefficients(z_aug, 1:this.s + 2, [], z + 1));
        end
        
        function name = get.name(this)
            name = ['FIRK-Gauss', num2str(this.s)];
        end
        
        function description = get.description(this)
            description = sprintf('%ith order, %i stage stage, fully-implict Gauss method', 2*this.s, this.s);
        end
        
        function order = get.order(this)
            order = 2 * this.s;
        end
        
    end
    
    methods (Access = private)
    
        function ObservablesHaveChanged(this, varargin)
            this.non_zero_output_indices = find(this.b);
            % -- set up stat objects -----------------------------------------------------------------------------------
            num_stages = size(this.A, 1);
            this.nonlinear_solver.stats.reset(num_stages);
            this.nonlinear_solver.linear_solver.stats.reset(num_stages);
        end
        
    end
    
end