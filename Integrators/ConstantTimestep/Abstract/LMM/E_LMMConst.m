
% ======================================================================================================================
%  Abstract Class for Explicit Linear Multistep Methods: 
%
%       y_{n+1} = \sum_{j=1}^p a(j) y_{n-p+j} + h * \sum_{j=1}^{p} b(j) * f_{n-p+j}
%
%   Properties
%       a (vector) - solution coefficients:   a(j) * y_{n-p+j}   for j = 1 ... p
%       b (vector) - derivative coefficients: h b(j) * f_{n-p+j} for j = 1 ... p
%
%       Examples:
%           Explicit Euler  ->  a = [1],     b = [1 0]
%           AB2             ->  a = [0 1],   b = [-1/2 3/2 0]
% ======================================================================================================================

classdef E_LMMConst < DI_LMMConst
    
    methods
        
        function this = E_LMMConst(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@DI_LMMConst(options);
            this.verifyCoefficientMatrices();
        end
        
    end
    
    methods (Access = protected)
        
        function verifyCoefficientMatrices(this)
            verifyCoefficientMatrices@DI_LMMConst(this);
            if(this.b(end) ~= 0)
                error('final coefficient of vector "b" must be of equal to zero for an explicit LMM.');
            end
        end
                
    end
    
end
