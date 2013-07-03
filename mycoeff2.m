function result = mycoeff2( exp , term)
%MYCOEFF Summary of this function goes here
%   Detailed explanation goes here
c = coeffs(expand(exp)+0.01 + term + term.^2,term);
if size(c,2) >= 3
    result = c(3) - 1;
else
    result = sym(0);
end

end

