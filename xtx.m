function ret = xtx(X)
% Convenience function for calculating the matrix product X'X.
%
% Syntax:
%   ret = xtx(X)
%
% Inputs:
%   X - [m,n] size, [double] type, input matrix
%
% Outputs:
%   ret - [n,n] size, [double] type, output matrix equal to X'X

ret = X'*X;
end
