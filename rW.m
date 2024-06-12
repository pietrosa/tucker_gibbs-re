function W = rW(kappa, m)
% A (modified) implementation of the rW function in the rstiefel package for R,
%   originally developed by Hoff. Used to support the rvmf function.
% The following description is modified from the rstiefel documentation.
%
% Auxilliary variable simulation for rejection sampling of rmf.vector, as described
%   in Wood (1994).
%
% Syntax:
%   W = rW(kappa, m)
%
% Inputs:
%   kappa - [1] size, [double] type, a positve scalar
%   a - [1] size, [int] type, a positive integer
%
% Outputs:
%   ret - [1] size, [double] type, scalar in (0,1) draws from the "W" distribution
%       described in Wood (1994)

b = (-2*kappa + sqrt(4*kappa^2 + (m-1)^2)) / (m-1);
x0 = (1-b)/(1+b);
c = kappa*x0 + (m-1)*log(1-x0^2);
accept = 0;
while accept==0
    Z = betarnd(0.5*(m-1), 0.5*(m-1));
    W = (1-(1+b)*Z) / (1-(1-b)*Z);
    U = rand();
    if kappa*W+(m-1)*log(1-x0*W)-c > log(U)
        accept = 1;
    end
    end
end
