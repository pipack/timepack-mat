
classdef ExponentialIntegratorConst < handle
    
    properties
        phi_evaluator
    end
    
    methods
        
        function this = ExponentialIntegratorConst(options)
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'phi_evaluator', KIOPS()}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            this.phi_evaluator = options.phi_evaluator;
        end
        
        function set.phi_evaluator(this, val)
            if(isa(val, 'PhiEvaluator'))
                this.phi_evaluator = val;
            else
                error('phi_evaluator must be of type PhiEvaluator');
            end
        end
        
    end
    
end
