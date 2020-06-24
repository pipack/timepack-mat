
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

classdef E_BLKConst < DI_BLKConst
    
    methods
        
        function this = E_BLKConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@DI_BLKConst(options);
        end
        
    end
    
    methods (Access = protected)
        
        function verifyCefficientMatrices(this)
            isstril = @(A) istril(A) && ~any(diag(A)); % tests if A is strictly lower triangular
            if(~isstril(this.D))
                error('Explicit Block Method must have strictly lower triangular D Matrix');
            end
            if(~isstril(this.C))
                error('Explicit Block Method must have strictly lower triangular C Matrix');
            end
        end
        
    end
    
end