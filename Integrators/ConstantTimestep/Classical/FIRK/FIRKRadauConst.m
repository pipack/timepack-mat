classdef FIRKRadauConst < FI_RKConst
    
    properties
        graphLineStyle = {};
        eval_RHS = true;
        extrapolate_initial_guess = false;
    end
    
    properties(SetObservable)
        s = 2; % number of stages (user can dyamically set)
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
        
        function this = FIRKRadauConst(options)
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
            z = this.nodes(this.s);
            ab = [zeros(this.s, 1), z(:)];
            A = double(iPolyCoefficients(z, ab, 'vpa'));
        end
        
        function b = get.b(this)
            z = this.nodes(this.s);
            ab = [0, 1];
            b = double(iPolyCoefficients(z, ab, 'vpa'));
        end
        
        function c = get.c(this)
            c = this.nodes(this.s);
        end
        
        function Ae = get.A_extrapolate(this)
            z = this.nodes(this.s);
            z_aug = [0; z]; % input, stage
            Ae = double(polyCoefficients(z_aug, 1:this.s + 1, [], z + 1));
            Ae(:,end+1) = 0; % add zero coefficinets for output (since output and final stage are equivalent)
        end
        
        function name = get.name(this)
            name = ['FIRK-Radau', num2str(this.s)];
        end
        
        function description = get.description(this)
            description = sprintf('%ith order, %i stage stage, fully-implict Radau method', 2*this.s - 1, this.s);
        end
        
        function order = get.order(this)
            order = 2 * this.s - 1;
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
        
        function n = nodes(this, s)
            if(s > 7) % use chebfun if q > 7
                if(isempty(which('radaupts')))
                    error('s > 7 requires Chebfun (https://www.chebfun.org/)');
                end
                n = (-flip(radaupts(s)) + 1) / 2;
                return;
            end
            
            n_all = {
                [1]'
                [0.333333333333333,   1]'
                [0.155051025721682,   0.644948974278318,  1]'
                [0.088587959512704,   0.409466864440735,  0.787659461760847,  1]'
                [0.057104196114518,   0.276843013638124,  0.583590432368917,  0.860240135656219,  1]'
                [0.039809857051469,   0.198013417873608,  0.437974810247386,  0.695464273353636,  0.901464914201174,  1]'
                [0.029316427159785,   0.148078599668484,  0.336984690281154,  0.558671518771550,  0.769233862030055,  0.926945671319741,  1]'
            };
            n = n_all{s};
        end
        
    end
    
end