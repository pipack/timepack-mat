classdef RKConst < IntegratorConst
    
    properties(SetAccess = protected)
        starting_times = 0;
    end
    
    methods
        
        function this = RKConst(options)
            this = this@IntegratorConst(options);
        end
        
    end
    
    methods (Access = protected)
        
        function indices = nonzeroStageIndices(this, A)
            
            num_stages = size(A, 1);
            indices    = cell(num_stages,1);
            for i = 1 : num_stages
                indices{i} = find(A(i, 1:(i-1)) ~= 0);
            end
            
        end
        
        function indices = nearestStageIndices(this, c)
            
            num_stages = length(c);
            indices    = cell(num_stages,1);
            c_extd     = [0; c(:)]; % include zero for methods that may not have Y_1 = y_n
            for i = 1 : num_stages
                [~, nearest_index] = min(abs(c(i) - c_extd(1:i)));
                indices{i} = nearest_index - 1;
            end
            
        end

    end
    
end