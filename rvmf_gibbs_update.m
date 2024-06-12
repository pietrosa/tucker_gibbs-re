function ret = rvmf_gibbs_update(M, X)
% A (simplified) implementation of rmf.matrix.gibbs in the rstiefel package for R,
%   originally developed by Hoff. In particular, this implementation forces
%   rscol=1 (ie, columns are updated one at a time).
%   See https://cran.r-project.org/web/packages/rstiefel/index.html
%   The following description is a modification of the rstiefel documentation.
%
% Simulate a random semi-orthogonal matrix from the matrix von Mises-Fisher
%   distribution using Gibbs sampling. This function only performs one Gibbs scan.
%
% Syntax:
%   ret = rmvf_gibbs_update(M, X)
%
% Inputs:
%   M - [n,p] size, [double] type, matrix parameter for the von Mises-Fisher distribution
%   X - [n,p] size, [double] type, current value of the random semi-orthogonal matrix
%
% Outputs:
%   ret - [n,p] size, [double] type, a new value of the input X obtained by Gibbs sampling

[U,S,V] = svd(M, 'econ');
H = U*S;
Y = X * V;
m = size(H,1);
R = size(H,2);
for iter = 1:R
    r = randsample(R,1);
    N = nullC(Y(:,setdiff(1:R,r)));
    y = rvmf(N'*H(:,r));
    Y(:,r) = N*y;
end
ret = Y*V';
end
