classdef Phipm < PhiEvaluator
    %KIOPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        stats = PhiStats() % IMPLEMENT STATS
        tolerance;
        krylov_max_dimension;
        krylov_min_dimension;
        krylov_start_dimension;
        symmetric;
        
    end
    
    methods
        
        function this = Phipm(options)
            %KIOPS Construct an instance of this class
            %   Detailed explanation goes here
            if(nargin == 0)
                options = struct();
            end
            options = this.DefaultOptions(options);
            this.tolerance = options.tolerance;
            this.krylov_max_dimension = options.krylov_max_dimension;
            this.krylov_min_dimension = options.krylov_min_dimension;
            this.krylov_start_dimension = options.krylov_start_dimension;
        end
        
        function [y, exit_flag] = phi(this, n, t, A)
            error('KIOPS does not support phi');
        end
        
        function [y, exit_flag] = phis(this, n, t, A)
            error('KIOPS does not support phis');
        end
        
        function [y, exit_flag] = phib(this, n, t, A, b)
            % Initialize phi_n(t*A) b
            y = this.phipm_raw(t, A, [zeros(size(b,1), n) b], this.tolerance, this.krylov_start_dimension, this.krylov_min_dimension, this.krylov_max_dimension, this.symmetric, true);
            exit_flag = true; % add functionality
        end
        
        function [y, exit_flag] = phibs(this, t, A, bs)
            % Initialize phi_0(t*A) b(:, 1) + phi_1(t * A) b(:, 2) + ... + phi_{length(b)-1}(t * A) b(:,end)
            y = this.phipm_raw(t, A, bs, this.tolerance, this.krylov_start_dimension, this.krylov_min_dimension, this.krylov_max_dimension, this.symmetric, true);
            exit_flag = true; % add functionality
        end
        
        function [y, exit_flag] = solve(this, t, A, bs)
            % compute solution to the differential equation
            %       y' = A y + bs(:, 1) + t * bs(:, 2) / 1! + t^2 * bs(:, 3) / 2! + ...
            % by initialize phi_0(t*A) b(:, 1) + t phi_1(t * A) b(:, 2) + ... + tphi_{length(b)-1}(t * A) b(:,end)
            y = this.phipm_raw(t, A, bs, this.tolerance, this.krylov_start_dimension, this.krylov_min_dimension, this.krylov_max_dimension, this.symmetric, false);
            exit_flag = true; % add functionality
        end
        
    end
    
    methods(Access = private)
        
        function options = DefaultOptions(this, options)
            if(~isfield(options, 'tolerance'))
                options.tolerance = 1e-10;
            end
            if(~isfield(options, 'krylov_max_dimension'))
                options.krylov_max_dimension = 128;
            end
            if(~isfield(options, 'krylov_min_dimension'))
                options.krylov_min_dimension = 8;
            end
            if(~isfield(options, 'krylov_start_dimension'))
                options.krylov_start_dimension = 8;
            end
        end
        
        function [w, stats] = phipm_raw(this, t, A, u, tol, m_init, mmin, mmax, symm, trimcoef)
            % PHIPM - Evaluates a linear combinaton of the phi functions
            %         evaluated at tA acting on vectors from u, that is
            %
            %         w = phi_0(tA) u(:, 1) + t phi_1(tA) u(:, 2) +
            %             t^2 phi_2(tA) u(:, 3) + ...
            %
            %         The evaluation expresses eveything in terms of the highest
            %         order phi function and evaluates the action of this on a
            %         vector using a Krylov technique and then computes w using
            %         the recurrence relation.
            %
            %         The size of the Krylov subspace is changed dynamically
            %         during the integration. The Krylov subspace is computed
            %         using Arnoldi if A is non-symmetric and Lancozs if A is
            %         symmetric.
            %
            %
            %
            % PARAMETERS:
            %   t    - constant value represent time.
            %   A    - the matrix argument of the phi functions.
            %   u    - the matrix with columns representing the vectors to be
            %          multiplied by the phi functions.
            %   tol  - the convergence tolarance required.
            %   symm - true if the matrix A is symmetric.
            %   m    - an estimate of the appropriate Krylov size.
            %
            % RETURNS:
            %   w        - the linear combination of the phi functions
            %              evaluated at tA acting on the vectors from u.
            %   stats(1) - number of substeps
            %   stats(2) - number of rejected steps
            %   stats(3) - number of Krylov steps
            %   stats(4) - number of matrix exponentials
            
            
            persistent V int nnze
            V = [];
            int = 0;
            
            % n is the size of the original problem
            % p is the number of phi functions
            [n, p] = size(u);
            Aisahandle = isa(A, 'function_handle');
            if ~Aisahandle
                nnze = nnz(A);
            else
                nnze = 10 * n; % wild guess
            end
            
            % Add extra column of zeros if p=1
            if p == 1
                p = 2;
                u = [u, zeros(size(u))];
            end
            
                        % Check inputs
            if nargin < 8
                trimcoef = false;
                if nargin < 7
                    mmax = 128;
                    if nargin < 6
                        mmin = 10;
                        if nargin < 5
                            m_init = mmin;
                            if nargin < 4
                                tol = 1.0e-7;
                                if nargin < 3
                                    error('Not enough input arguments.');
                                end
                            end
                        end
                    end
                end
            end
            m_new = m_init;
            
            % Preallocate matrices
            % % V = zeros(n, mmax + 1);
            % % int = zeros(n, p);
            % if isempty(V) && isempty(int)
            %   V = zeros(n, mmax + 1);
            %   int = zeros(n, p);
            % elseif numel(int) ~= n * p
            %   int = zeros(n, p);
            % elseif numel(V) ~= n * (mmax + 1)
            %   V = zeros(n, mmax + 1);
            %   int = zeros(n, p);
            % end
            if isempty(V) && isempty(int)
                V = zeros(n, mmax+1);
                int = zeros(n, p);
            elseif numel(V) ~= n*(mmax+1)
                V = zeros(n, mmax+1);
                int = zeros(n, p);
            elseif numel(int) ~= n*p
                int = zeros(n, p);
            end
            
            if isempty(nnze)
                nnze = nnz(A);
            end
            
            % Initializing the variables
            step = 0;
            krystep = 0;
            ireject = 0;
            reject = 0;
            exps = 0;
            happy = 0;
            sgn = sign(t(end));
            t_now = 0;
            t_out = abs(t(end));
            j = 0;
            sizet = size(t);
            numSteps = sizet(2);
            
            % Compute and initial starting approximation for the timestep
            tau = t_out;
            
            % Setting the safety factors and tolarance requirements
            gamma = 0.8;
            delta = 1.2;
            
            % Used for toeplitz trick
            cidx = (0:p-1)';
            ridx = p:-1:1;
            idx = cidx(:,ones(p,1)) + ridx(ones(p,1),:);
            
            % Initial condition
            w = zeros(length(u(:,1)), numSteps);
            w(:,1) = u(:, 1);
            oldm = NaN; oldtau = NaN; omega = NaN;
            orderold = true; kestold = true;
            
            l = 1;
            while t_now < t_out
                
                % Compute necessary starting information
                if j == 0
                    
                    % Initialize the matrices V and H
                    H = zeros(mmax + p, mmax + p);
                    x = [zeros(1, p - 1), cumprod([1, t_now ./ (1:p - 1)])];
                    up = u * x(idx);  % Code inspired by toeplitz.m
                    %up2 = u * x(idx);
                    
                    % Compute the update factors
                    int(:, 1) = w(:,l);
                    for i = 1:p - 1
                        
                        % Determine if matrix free
                        if ~Aisahandle
                            int(:, i + 1) = A * int(:, i) + up(:, i + 1);
                        else
                            int(:, i + 1) = A(int(:, i)) + up(:, i + 1);
                        end
                        
                    end
                    
                    % Normalize initial vector
                    beta = norm(int(:, end));
                    if beta == 0
                        
                        % Multiplying with a zero vector, hence result is zero
                        % Finish all in one step
                        reject = reject + ireject;
                        step = step + 1;
                        tau = t_out - t_now;
                        w(:,l) = w(:,l) + int(:, 2:p - 1) * cumprod(tau * 1 ./ (1: p - 2)');
                        break;
                        
                    end
                    
                    % The first Krylov basis vector
                    V(:, 1) = int(:, end) ./ beta;
                    
                end
                
                % Check if matrix is symmetric
                if symm
                    
                    % Symmetric use the Lanczos process
                    while j < m_init
                        
                        % Determine if matrix free
                        j = j + 1;
                        if ~Aisahandle
                            vv = A * V(:, j);
                        else
                            vv = A(V(:, j));
                        end
                        H(j, j) = V(:, j)' * vv;
                        if j == 1
                            vv = vv - H(j, j) * V(:, j);
                        else
                            vv = vv - H(j - 1, j) * V(:, j - 1) - H(j, j) * V(:, j);
                        end
                        krystep = krystep + 1;
                        s = norm(vv);
                        
                        % Happy breakdown
                        if s < tol
                            happy = 1;
                            tau = t_out - t_now;
                            break;
                        end
                        H(j + 1, j) = s;
                        H(j, j + 1) = s;
                        V(:, j + 1) = vv ./ s;
                        
                    end
                    
                    % Keep a record of H
                    H2 = H;
                    H(m_init, m_init + 1) = 0;
                    
                    % Matrix is not symmetric
                else
                    
                    % Not symmetric use the Arnoldi process
                    while j < m_init
                        
                        % Determine if matrix free
                        j = j + 1;
                        if ~Aisahandle
                            vv = A * V(:, j);
                        else
                            vv = A(V(:, j));
                        end
                        for i = 1:j
                            H(i, j) = V(:, i)' * vv;
                            vv = vv - H(i, j) * V(:, i);
                        end
                        krystep = krystep + 1;
                        s = norm(vv);
                        
                        % Happy breakdown
                        if s < tol
                            happy = 1;
                            tau = t_out - t_now;
                            break;
                        end
                        H(j + 1, j) = s;
                        V(:, j + 1) = vv ./ s;
                        
                    end
                    
                    % Keep a record of H
                    H2 = H;
                    
                end
                
                % We use the vector e1 in the computations
                H(1, j + 1) = 1;
                
                % Construct the augmented matrix
                for i = 1:p - 1
                    H(j + i, j + i + 1) = 1;
                end
                h = H(j + 1, j);
                H(j + 1, j) = 0;
                
                % Compute the exponential of the augmented matrix
                
                [F, hnorm] = this.expm_raw(sgn * tau * H(1:j + p, 1:j + p));
                exps = exps + 1;
                % Local truncation error estimation
                err = abs(beta * h * F(j, j + p));
                %err = norm(beta * h * F(j, j + p - 1) * V(:, j+1));
                
                % Error per unit step
                oldomega = omega;
                omega = t_out * err / (tau * tol);
                
                % Estimate order
                if m_init == oldm && tau ~= oldtau && ireject >= 1
                    order = max(1, log(omega/oldomega) / log(tau/oldtau));
                    orderold = false;
                elseif orderold || ireject == 0
                    orderold = true;
                    order = j/4;
                else
                    orderold = true;
                end
                % Estimate k
                if m_init ~= oldm && tau == oldtau && ireject >= 1
                    kest = max(1.1, (omega/oldomega) ^ (1/(oldm-m_init)));
                    kestold = false;
                elseif kestold || ireject == 0
                    kestold = true;
                    kest = 2;
                else
                    kestold = true;
                end
                
                % This if statement is the main difference between fixed and
                % variable m
                oldtau = tau; oldm = m_init;
                if happy == 1
                    
                    % Happy breakdown; wrap up
                    omega = 0;
                    tau_new = tau;
                    m_new = m_init;
                    
                elseif j == mmax && omega > delta
                    
                    % Krylov subspace to small and stepsize to large
                    tau_new = tau * (omega / gamma) ^ (-1 / order);
                    
                else
                    
                    % Determine optimal tau and m
                    tau_opt = tau * (omega / gamma) ^ (-1 / order);
                    m_opt = max(1, ceil(j + log(omega / gamma) / log(kest)));
                    nom = 5 + max(log(hnorm), 0) / log(2); % number of mult's in expm
                    
                    if symm
                        
                        % Cost of Lanczos; a factor of 2 has been ignored
                        cost1 = ((j + p) * nnze + 3 * (j + p) * n ...
                            + nom * (j + p - 1)^3) * ceil((t_out-t_now) / tau_opt);
                        cost2 = ((m_opt + p) * nnze + 3*(m_opt + p) * n ...
                            + nom * (m_opt + p - 1)^3) * ceil((t_out-t_now) / tau);
                    else
                        
                        % Cost of Arnoldi
                        cost1 = ((j + p) * nnze + (j^2 + 3 * p + 2) * n ...
                            + nom * (j + p - 1)^3) * ceil((t_out-t_now) / tau_opt);
                        cost2 = ((m_opt + p) * nnze + (m_opt^2 + 3 * p + 2) * n ...
                            + nom * (m_opt + p - 1)^3) * ceil((t_out-t_now) / tau);
                        
                    end
                    
                    % Determine whether to vary tau or m
                    if cost1 < cost2
                        tau_new = tau_opt;
                        m_new = m_init;
                    else
                        m_new = m_opt;
                        tau_new = tau;
                    end
                    
                end
                
                % Check error against target
                if omega <= delta
                    
                    blownTs = 0;
                    nextT = t_now + tau;
                    for k = l:numSteps
                        if t(k) < nextT
                            blownTs = blownTs + 1;
                        end
                    end
                    
                    % Copy current w to w we continue with.
                    w(:,l + blownTs) = w(:,l);
                    
                    for k = 0:blownTs - 1
                        tauPhantom = t(l+k) - t_now;
                        F2 = this.expm_raw(sgn * tauPhantom * H(1:j + p, 1:j + p));
                        up2 = w(:,l+blownTs) + int(:, 2:p - 1) * cumprod(tauPhantom * 1 ./ (1: p - 2)');
                        F2(j + 1, j + p - 1) = h * F2(j, j + p);
                        w(:,l+k) = beta * V(:, 1:j + 1) * F2(1:j + 1, j + p - 1) + up2;
                    end
                    
                    % Advance l.
                    l = l + blownTs;
                    
                    % Yep, got the required tolerance; update
                    reject = reject + ireject;
                    step = step + 1;
                    up = w(:,l) + int(:, 2:p - 1) * cumprod(tau * 1 ./ (1: p - 2)');
                    
                    % Using the corrected quantity
                    F(j + 1, j + p - 1) = h * F(j, j + p);
                    w(:,l) = beta * V(:, 1:j + 1) * F(1:j + 1, j + p - 1) + up;
                    
                    % Update t
                    %fprintf('Step = %d, basisVectors = %d\n', step, j);
                    t_now = t_now + tau;
                    j = 0;
                    ireject = 0;
                    
                else
                    
                    % Nope, try again
                    H = H2;
                    ireject = ireject + 1;
                    
                end
                
                % Safety factors for tau
                %tau = min(t_out - t_now, max(tau/5, min(5*tau, tau_new)));
                tau = min(t_out - t_now, max(tau/5, min(2*tau, tau_new)));
                
                % Safety factors for m
                m_init = max(1, min(mmax, max(floor(3/4*m_init), min(m_new, ceil(4/3*m_init)))));
                
            end
            
            if t(1)~=1 && trimcoef
                if length(t)==1
                    w(:,l) = w(:,l)*(1/t(l))^p;
                else
                    phiHan = find(max(abs(u)));
                    
                    if isempty(phiHan)
                        phiHan = length(u);
                    else
                        phiHan = phiHan-1;
                    end
                    
                    for l = 1:numSteps
                        %w(:,l) = w(:,l)*(1/t(l).^phiHan);
                    end
                end
            end
            
            
            stats = [step, reject, krystep, exps];
            
        end
        
        function [F, normA] = expm_raw(this, A)
            %   Compute the matrix exponential of A using a scaling and squaring
            %   algorithm with a Pade approximation.
            %
            %   Reference:
            %   N. J. Higham, The scaling and squaring method for the matrix
            %   exponential revisited. SIAM J. Matrix Anal. Appl.,
            %   26(4) (2005), pp. 1179-1193.
            %
            switch class(A)
                case 'double'
                    m_vals = [3 5 7 9 13];
                    % theta_m for m=1:13.
                    theta = [%3.650024139523051e-008
                        %5.317232856892575e-004
                        1.495585217958292e-002  % m_vals = 3
                        %8.536352760102745e-002
                        2.539398330063230e-001  % m_vals = 5
                        %5.414660951208968e-001
                        9.504178996162932e-001  % m_vals = 7
                        %1.473163964234804e+000
                        2.097847961257068e+000  % m_vals = 9
                        %2.811644121620263e+000
                        %3.602330066265032e+000
                        %4.458935413036850e+000
                        5.371920351148152e+000];% m_vals = 13
                case 'single'
                    m_vals = [3 5 7];
                    % theta_m for m=1:7.
                    theta = [%8.457278879935396e-004
                        %8.093024012430565e-002
                        4.258730016922831e-001  % m_vals = 3
                        %1.049003250386875e+000
                        1.880152677804762e+000  % m_vals = 5
                        %2.854332750593825e+000
                        3.925724783138660e+000];% m_vals = 7
                otherwise
                    error(message('MATLAB:expm:inputType'))
            end
            
            normA = norm(A,1);
            if normA <= theta(end)
                
                % no scaling and squaring is required.
                for i = 1:length(m_vals)
                    if normA <= theta(i)
                        F = this.epade_raw(m_vals(i),A);
                        break;
                    end
                end
                
            else
                
                % Scale the matrix A so that a Padé approximant of degree 13 will be accurate
                [t, s] = log2(normA/theta(end));
                s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
                
                A  = A*(1/2^s);    % Scaling
                
                F = this.epade_raw(m_vals(end),A);
                
                for i = 1:s
                    F = F*F;  % Squaring
                end
                
            end
            
        end
        
        function F = epade_raw(~, m,A)
            %   Pade' approximant to exponential.
            %   Compute the degree M diagonal Pade' approximant to EXP(A),
            %   where M = 3, 5, 7, 9 or 13.
            %   Series are evaluated in decreasing order of powers, which is
            %   in approx. increasing order of maximum norms of the terms.
            n = length(A);
            I = eye(n,class(A));
            A2 = A*A;
            % Evaluate Pade approximant.
            switch m
                case 3
                    c = [120, 60, 12, 1];
                    U = A * (c(4)*A2 + c(2)*I);
                    V = c(3)*A2 + c(1)*I;
                case 5
                    c = [30240, 15120, 3360, 420, 30, 1];
                    A4 = A2*A2;
                    U = A * (c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = c(5)*A4 + c(3)*A2 + c(1)*I;
                case 7
                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
                    A4 = A2*A2;
                    A6 = A2*A4;
                    U = A * (c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I;
                case 9
                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
                        2162160, 110880, 3960, 90, 1];
                    A4 = A2*A2;
                    A6 = A2*A4;
                    A8 = A4*A4;
                    U = A * (c(10)*A8 + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = c(9)*A8 + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I;
                case 13
                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
                        1187353796428800,  129060195264000,   10559470521600,   ...
                        670442572800,      33522128640,       1323241920,       ...
                        40840800,          960960,            16380,  182,  1];
                    A4 = A2*A2;
                    A6 = A2*A4;
                    U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*I);
                    V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*I;
            end
            F = (V-U)\(2*U) + I;
        end
        
    end
    
end