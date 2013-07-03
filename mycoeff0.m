function result = mycoeff0( exp , term)
%MYCOEFF Summary of this function goes here
%   Detailed explanation goes here
c = coeffs(expand(exp) + 0.01,term);
if size(c,2) >= 1
    result = c(1) -  0.01;
else
    result = sym(0);
end

end

