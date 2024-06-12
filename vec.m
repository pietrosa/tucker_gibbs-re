function v = vec(x)
% Vectorizes a given matrix or tensor
%
% Syntax:
%   v = vec(x) is the vectorization of x
%
% Inputs:
%   x - [n1,...,nq] size (q arbitrary), [double|tensor] type, matrix/tensor to
%       be vectorized
%
% Outputs:
%   v - [n1*...*nq] size, [double] type, vectorization of x

v = reshape( x, numel( x ), 1 );
end
