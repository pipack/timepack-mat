
% ======================================================================================================================
%  Explicit Block Method
%       y^{(n+1)} = A * y^{(n)} + r * B * y^{(n+1)} + C * y^{(n+1)} + r * D * f^{[n+1]}
%  with optional output point
%       y_out = a_out * y^{(n)} + r * b_out * y^{(n+1)} + c_out * y^{(n+1)} + r * d_out * f^{[n+1]} + r * e_out * f_out
%
%   A parallel diagonally implicit block method must satisfy the following conditions:
%       1. The matrix D must be a the zero matrix.
%       2. The matrix C must be strictly lower triangular
% ======================================================================================================================

classdef E_PBLKConst < DI_BLKConst
    
    methods
        
        function this = E_PBLKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@DI_BLKConst(options);
        end
        
    end
    
    methods (Access = protected)
        
        function verifyCefficientMatrices(this)
            if(any(this.D))
                error('Coefficient Matrix D must be empty for an explicit parallel block method.');
            end
            if(any(this.C))
                error('Coefficient Matrix C must be empty for an explicit parallel block method.');
            end
        end
        
    end
    
end