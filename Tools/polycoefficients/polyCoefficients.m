function [C, FD] = polyCoefficients(z, A, B, tau)
%COEFFICIENT Summary of this function goes here
%   A locations where value matches
%   B locations where derivative matches

len_A = length(A);
len_B = length(B);
p     = len_A + len_B;

M = sym(zeros(p));
% -- Value Conditions --------
for i=1:len_A
    z0 = z(A(i));
    for j=1:p
        M(i,j) = z0^(j-1);
    end
end
% -- Derivative Conditions ---
for i=1:len_B
    z0 = z(B(i));
    for j=2:p
        M(i+len_A,j) = (j-1) * z0^(j-2);
    end
end
% -- Evaluation --------------
if(nargin ==  4)
    E = sym(zeros(length(tau), p));
    for i=1:length(tau)
        for j=1:p
            E(i,j) = tau(i)^(j-1);
        end
    end
    C = transpose(transpose(M) \ transpose(E)); % equivalent and more stable than C = E*inv(M);
else
    C = inv(M);
end
if(nargout == 2)
    FD = inv(M);
end
end