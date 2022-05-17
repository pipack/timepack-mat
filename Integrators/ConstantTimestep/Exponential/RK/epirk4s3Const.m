classdef epirk4s3Const < IntegratorConst & ExponentialIntegratorConst
    properties
        graph_line_style = {};
        matrix_free = false;
    end
    
    properties (SetAccess = protected)
        order = 4;
        name = 'epirk4s3';
        description = 'epirk4s3'
        starting_times = 0;
        
    end
        
    methods
        function this = epirk4s3Const(options)
            if(nargin == 0)
                options = struct();
            end
            this = this@IntegratorConst(options);
            this = this@ExponentialIntegratorConst(options);
        end
        
    end
    
    methods (Access = protected)
        
        function [step_struct] = initStepStruct(this, t_in, y_in, problem)
            step_struct = struct();
        end
        
        function [t_out, y_out, step_struct] = step(this, t_in, y_in, step_struct, problem)
            
            %  Coefficients for
            %  U_2     = u_n + alpha21*(p211 Phi_1(g21 z)+p212 Phi_2(g21 z)+p213 Phi_3(g21 z))hF
            %  U_3     = u_n + alpha31*(p311 Phi_1(g31 z)+p312 Phi_2(g31 z)+p313 Phi_3(g31 z))hF + alpha32*Psi(g32 z) hR(U2)
            %  u_{n+1} = u_n + Phi_1(z)hF + (b2p3 Phi_3(z) + b2p4 Phi_4(z)) hR(U2) + (b3p3 Phi_3(z) + b3p4 Phi_4(z)) hR(U3)
            
            mtrx_free = this.matrix_free;
            
            g21 = 1 / 8;
            g31 = 1 / 9;
            
            alpha21 = 1 / 8;
            alpha31 = 1 / 9;
            
            b2p3 = -1024;
            b2p4 = 27648;
            
            b3p3 = 1458;
            b3p4 = -34992;
            
            % Setup initial state.
            h = this.h;
            
            % Correctly order the g-coefficients for adaptive-Krylov
            N = length(y_in);
            zeroVec = zeros(N,1);
            gCoeffVec = [g31 g21];
            KryIndex = [2 1];
            
            step_start_time = tic;
            hF = h * problem.RHS(y_in);  % calculate right hand side
            
            if(mtrx_free)
                hA = @(x) h * problem.Jx(y_in, x);  % function handle for product
            else
                hA = h * problem.J(y_in);  % calculate full Jacobian
            end
            
            % -- Stage 1 -------------------------------------------------------------------------------------------
            [temp1, clean_exit] = this.phi_evaluator.phibs(gCoeffVec, hA, [zeroVec, hF]);
            if(~clean_exit)
                t_out = NaN;
                y_out = NaN;
                return;
            end
            U2 = y_in + alpha21*temp1(:,KryIndex(1));
            if(mtrx_free)
                hb1 = h * problem.RHS(U2) - hF - hA(U2 - y_in);     % Calculate residual r(U2)
            else
                hb1 = h * problem.RHS(U2) - hF - hA*(U2 - y_in);     % Calculate residual r(U2)
            end
            
            % -- Stage 2 -------------------------------------------------------------------------------------------
            U3 = y_in + alpha31*temp1(:,KryIndex(2));
            if(mtrx_free)
                hb2 = h * problem.RHS(U3) - hF - hA(U3 - y_in);    % Calculate residual r(U3)
            else
                hb2 = h * problem.RHS(U3) - hF - hA*(U3 - y_in);    % Calculate residual r(U3)
            end
            
            % -- Stage 3 -------------------------------------------------------------------------------------------
            [temp3, clean_exit] = this.phi_evaluator.phibs(1, hA, [zeroVec, hF, zeroVec, b2p3*hb1+b3p3*hb2, b2p4*hb1+b3p4*hb2]);
            if(~clean_exit)
                t_out = NaN;
                y_out = NaN;
                return;
            end
            
            y = y_in + temp3(:,1);
            
            % Advance simulation state.
            t_out = t_in + h;
            y_out = y;
            
            % Update counters.
            this.step_stats.recordStep(toc(step_start_time));
            
        end
        
        function [t_user, y_user] = userOutput(this, t_in, y_in, struct_in, t_out, y_out, struct_out, problem)
            t_user = t_out;
            y_user = y_out;
        end
        
    end
    
end