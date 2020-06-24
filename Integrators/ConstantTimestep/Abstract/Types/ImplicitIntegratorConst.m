
classdef ImplicitIntegratorConst < handle
    
    properties
        nonlinear_solver
    end
    
    methods
        
        function this = ImplicitIntegratorConst(options)
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'nonlinear_solver', Newton()}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            this.nonlinear_solver = options.nonlinear_solver;
        end
        
        function set.nonlinear_solver(this, val)
            if(isa(val, 'NonlinearSolver'))
                this.nonlinear_solver = val;
            else
                error('phi_evaluator must be of type PhiEvaluator');
            end
        end
        
    end
    
end
