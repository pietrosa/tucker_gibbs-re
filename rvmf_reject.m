function ret = rvmf_reject(M)
% A (modified) implementation of the rejection sampler for the von Mises-Fisher
%   distribution, originally implemented as rmf.matrix (by Hoff) in the rstiefel
%   package for R.
%   The following description is modified from the rstiefel documentation.
%
% Simulate a random orthonormal matrix from the von Mises-Fisher distribution.
%
% Syntax:
%   ret = rvmf_reject(M)
%
% Inputs:
%   M - [n,p] size, [double] type, matrix argument of the von Mises-Fisher distribution
%
% Outputs:
%   ret - [n,p] size, [double] type, simulated draw from the von Mises-Fisher distribution

if size(M,2)==1
    ret = rvmf(M);
else
    [U,S,V] = svd(M);
    H = U*S;
    m = size(H,1);
    R = size(H,2);
    cmet = false;
    rej = 0;
    while ~cmet
        U = zeros(m,R);
        U(:,1) = rvmf(H(:,1));
        lr = 0;
        for j = 2:R
            N = nullC(U(:,1:(j-1)));
            x = rvmf(N'*H(:,j));
            U(:,j) = N*x;
            if S(j,j)>0
                xn = sqrt(sum((N'*H(:,j)).^2, 'all'));
                xd = sqrt(sum(H(:,j).^2));
                lbr = log(besseli(0.5*(m-j-1), xn, 1)) - log(besseli(0.5*(m-j-1), xd, 1));
                if isnan(lbr)
                    lbr = 0.5*(log(xd)-log(xn));
                end
                lr = lr + lbr + (xn-xd) + 0.5*(m-j-1)*(log(xd)-log(xn));
            end
        end
        cmet = log(rand()) < lr;
        rej = rej + (1-1*cmet);
    end
    ret = U*V';
end
