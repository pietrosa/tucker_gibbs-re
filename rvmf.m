function ret = rvmf(kv)
% A (modified) implementation of rmf.vector in the rstiefel package for R,
%   originally developed by Hoff.
%   The following description is modified from the rstiefel documentation.
%
% Simulate a random normal vector from the von Mises-Fisher distribution as
%   described in Wood (1994).
%
% Syntax:
%   ret = rvmf(kv)
%
% Inputs:
%   kv - [n] size, [double] type, vector parameter for the (vector)
%   von Mises-Fisher distribution.
%
% Outputs:
%   ret - [n] size, [double] type, vector draws from the (vector) von Mises-Fisher
%       distribution with matrix parameter 


kappa = sqrt(sum(kv.^2));
mu = kv./kappa;
p = length(kv);
if kappa==0
    ret = randn(p,1);
    ret = ret / sqrt(sum(ret.^2));
else
    if p==1
        ret = 2*binornd(1, 1/(1+exp(2*kv)))-1;
    else
        if kappa>1e4
            ret = mu;
        else
            W = rW(kappa, p);
            v = randn(p-1,1);
            v = v/sqrt(sum(v.^2));
            x = [sqrt(1-W^2).*v; W];
            ret = [nullC(mu), mu] * x;
        end
    end
end

end
