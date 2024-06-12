function ret = nullC(M)
% An implementation of NullC in the rstiefel package for R, originally developed
%   by Hoff. See https://cran.r-project.org/web/packages/rstiefel/index.html
%   The following description is a modification of the rstiefel documentation.
%
% Given a matrix ‘M’, find a matrix ‘N’ giving a basis for the null space. This
%   is a modified version of Null from the package MASS.
%
% Syntax:
%   N = NullC(M) gives an orthonormal matrix N such that N'M is a matrix of zeros.
%
% Inputs:
%   M - [n,p] size, [double] type, input matrix
%
% Outputs:
%   ret - [n, n-rank(M)] size, [double] type, orthonorml matrix
%       such that N'M is a matrix of zeros.


[Q,~]  = qr(M);
rk = rank(M);
ret = Q(:,(rk+1):size(Q,2));
end
