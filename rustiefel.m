function ret = rustiefel(rc)
% Draws a single value from the uniform distribution U(V) on the Stiefel manifold.
%   By a well-known result, the polar part (Q) or the polar (QR) decomposition
%       of a matrix with iid mean-zero, normal entries follows U(V).
%
% Syntax:
%   ret = rustiefel([n, p]) draws a random n*p semi-orthogonal matrix from the
%       Stiefel manifold
%
% Inputs:
%   rc - [2] size, [int] type, vector containing the number of rows and columns
%       of the matrix to be drawn
%
% Outputs:
%   ret - [n,p] size, [double] type, n*p semi-orthogonal matrix drawn from the
%       uniform distribution on the Stiefel manifold

[ret,~] = qr(randn(rc), 'econ');
end
