classdef NonlinearSolver < handle
    %NONLINEARSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract = true)
        stats;
    end
    
    methods(Abstract = true)
        [y, exit_flag, residual] = solve(this)   % Solve system F(x) = 0
        [y, exit_flag, residual] = solveBC(this) % Solve system x = b + c * F(x)     (c scalar)
        
        [y, exit_flag, residual] = solveBCKronF(this);
        % Solve the system x = b + kron(C, F) where C is a matrix. For example, for 3x3 C we have:
        %  | x_1 |   | C_11 F(x_1), C_12 F(x_2), C_13 F(x_3) |    
        %  | x_2 | - | C_21 F(x_1), C_22 F(x_2), C_23 F(x_3) | = b  
        %  | x_3 |   | C_31 F(x_1), C_32 F(x_2), C_33 F(x_3) |    
        
    end
    
end