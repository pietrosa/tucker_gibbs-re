function ret = kron_diag(a, b, varargin)
% Computes the Kronecker product of two (or more) diagonal matrices.
%   Calculates A x B x C x ..., where "x" denotes the Kronecker product and
%   A, B, C, ... are diagonal matrices. Only diagonal elements should be
%   provided as inputs.
%
% Syntax:
%   out = kron_diag(a,b) calculates the Kronecker product diag(a) x diag(b)
%   out = kron_diag(a,b,c,...) calculates the Kronecker product diag(a) x diag(b) x diag(c) x ...
%
% Inputs:
%   a - [na] size, [double] type, diagonal elements of the first diagonal matrix
%   b - [nb] size, [double] type, diagonal elements of the second diagonal matrix
%   varargin - additional arguments with the same type as a, b (sizes may differ);
%       the diagonal elements of the third, fourth, etc. matrices in the
%       Kronecker product
%
% Outputs:
%   ret - [na*nb*...] size, [double] type, diagonal elements of the Kronecker
%       product diag(a) x diag(b) x ...

ret = repmat(b,1,length(a)) .* repelem(a,length(b));
if nargin>2
    ret = kron_diag(ret, varargin{:});
end

end
