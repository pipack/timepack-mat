
classdef LinearlyImplicitIntegratorConst < handle
    
    properties
        linear_solver
    end
    
    methods
        
        function this = LinearlyImplicitIntegratorConst(options)
            % -- set local properties ----------------------------------------------------------------------------------
            default_field_value_pairs = {
                {'linear_solver', MLDV()}
            };
            options = setDefaultOptions(options, default_field_value_pairs);
            % -- set props ---------------------------------------------------------------------------------------------
            this.linear_solver = options.linear_solver;
        end
                
        function set.linear_solver(this, val)
            if(isa(val, 'LinearSolver'))
                this.linear_solver = val;
            else
                error('linear_solver must be of type LinearSolver');
            end
        end
        
    end
    
end
