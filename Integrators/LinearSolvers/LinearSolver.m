classdef LinearSolver < handle
    %LINEARSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract = true)
        stats;
    end
    
    properties(Abstract = true, SetAccess = protected)
        allows_matrix_free % (bool) if true, solve and solveBC must accept a function_handle that computes A*x, and a standard matrix
    end
    
    methods(Abstract = true)
        [y, flag, residual] = solve(this);       % Solve system A * x = b
        [y, flag, residual] = solveBC(this);     % Solve system x = b + c * A * x    (c scalar)
        
        [y, flag, residual] = solveBCKronA(this);
        % Solve the system x = b + kron(C, A) where C is a matrix. For example, for 3x3 C we have:
        %       | C_11 A, C_12 A, C_13 A |    
        % ( I - | C_21 A, C_22 A, C_23 A | ) x = b  
        %       | C_31 A, C_32 A, C_33 A |    
        
        [y, flag, residual] = solveBCKronAs(this);
        % Solve the system kron(C,I) * blkdiag(As{:}) where size(C,2) = length(As). For example, for 3x3 C we have:
        %       | C_11 A{1}, C_12 A{2}, C_13 A{3} |    
        % ( I - | C_21 A{1}, C_22 A{2}, C_23 A{3} | ) x = b  
        %       | C_31 A{1}, C_32 A{2}, C_33 A{3} |     
        
    end
    
end

