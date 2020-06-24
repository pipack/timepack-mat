classdef LinearSolver < handle
    %NONLINEARSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract = true)
        stats;
    end
    
    properties(Abstract = true, SetAccess = protected)
        allows_matrix_free % (bool) if true, solve and solveBC must accept a function_handle that computes A*x, and a standard matrix
    end
    
    methods(Abstract = true)
        [y, flag, residual] = solve(this);       % Solve system A * x = b
        [y, flag, residual] = solveBC(this);     % Solve system of form x = b + c * A * x 
    end
    
end

