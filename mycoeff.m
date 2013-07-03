function result = mycoeff( exp , term)
%MYCOEFF Summary of this function goes here
%   Detailed explanation goes here
c = coeffs(expand(exp)+ 0.01 + term,term);
if size(c,2) >= 2
    result = c(2) - 1;
else
    result = sym(0);
end

end

