classdef RK4Const < E_RKConst
    
    properties
        graphLineStyle = {'kd-'};
    end
    
    properties(SetAccess = protected)
        name  = 'RK4';
        description = 'Classical Fourth Order Runge Kutta';
        order = 4;
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
        b = [1/6 1/3 1/3 1/6];
    end
    
end