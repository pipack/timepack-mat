classdef E_RKConst < RKConst
    
    properties
        graph_line_style = {};
    end
    
    properties(Abstract = true, SetAccess = protected)
        A % RK Matrix (Must be strictly lower triangular)
        b % RK b vector
    end
    
    properties(SetAccess = protected)
        non_zero_stage_indices  = {};
        non_zero_output_indices = [];
    end
    
    methods
        
        function this = E_RKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@RKConst(options);
            this.non_zero_stage_indices = this.nonzeroStageIndices(this.A);
            this.non_zero_output_indices = find(this.b);
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct] = initStepStruct(this, ~, y_in, ~)
            
            num_stages = size(this.A, 1);
            ode_dim    = length(y_in);
            
            step_struct = struct(                   ...
                'AT', transpose(this.A),            ...
                'bT', transpose(this.b),            ...
                's',  num_stages,                   ...
                'F',  zeros(ode_dim, num_stages)    ... % Stage Derivative Vector
                );
            
        end
        
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem, ~)
            
            step_start_time = tic; % -- start step time clock ----------------------------------------------------------
            
            h = this.h; % stepsize
            num_stages = step_struct.s;
            
            for j = 1:num_stages
                nzi = this.non_zero_stage_indices{j};
                start_rhs_time = tic;
                step_struct.F(:,j) = problem.RHS(y_in + h * step_struct.F(:, nzi) * step_struct.AT(nzi, j));
                this.rhs_stats.recordRHSEval(toc(start_rhs_time));
            end
            
            nzo   = this.non_zero_output_indices;
            y_out = y_in + step_struct.F(:, nzo) * (h * step_struct.bT(nzo));
            t_out = t_in + h;
            
            this.step_stats.recordStep(toc(step_start_time)); % -- stop step time clock --------------------------------
            
        end
        
    end
    
end