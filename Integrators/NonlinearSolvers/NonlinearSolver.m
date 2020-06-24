classdef NonlinearSolver < handle
    %NONLINEARSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract = true)
        stats;
    end
    
    methods(Abstract = true)
        [y, exit_flag, residual] = solve(this)       % Solve system F(x) = 0
        [y, exit_flag, residual] = solveBC(this)     % Solve system of form x = b + c * F(x) 
    end
    
end

