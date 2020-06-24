function [C, FD] = iPolyCoefficients(z, b, numeric_type)
%COEFFICIENT Summary of this function goes here
%   A locations where derivative matches
%   B locations where value matches

q = length(z);
num_ints = size(b,1);

if(strcmp(numeric_type, 'sym'))
    I = diag(sym(1)./sym(1:q));
    M = sym(zeros(q));
    E = sym(zeros(num_ints, q));
elseif(strcmp(numeric_type, 'vpa'))
    I = diag(vpa(1)./vpa(1:q));
    M = vpa(zeros(q));
    E = vpa(zeros(num_ints, q));
else
    I = diag(1./(1:q));
    M = zeros(q);
    E = zeros(num_ints, q);
end

for i=1:q
    for j=1:q
        M(i,j) = z(i)^(j-1);
    end
end


% -- Evaluation --------------
if(nargin >=  2)
    for i=1:num_ints
        for j=1:q
            E(i,j) = b(i,2)^j - b(i,1)^j;
        end
    end
    C = transpose(transpose(M) \ transpose(E*I)); % equivalent and more stable than C = E*I*inv(M);
else
    C = inv(M);
end
if(nargout == 2)
    FD = inv(M);
end
end