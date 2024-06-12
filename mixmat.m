function ret = mixmat(n, meth)
% Generates mixing matrices. Supports tucker_gibbs.
%
% Syntax:
%   ret = mixmat(n,'dct') gives an n*n orthonormal matrix HD, where H is the n*n
%       discrete cosine transform matrix, and D is an n*n diagonal Rademacher
%       matrix
%
% Inputs:
%   n - [1] size, [int] type, size of the desired output matrix
%   meth - [1] size, [string] type, type of mixing matrix to generate. Currently
%       only 'dct' is supported; any other value will return the identity matrix.
%
% Outputs:
%   ret - [n, n] size, [double] type, orthonormal mixing matrix (if meth='dct',
%       this is a DCT matrix with columns randomly multiplied by -1 or +1; otherwise,
%       this is the identity matrix)

if strcmp(meth,'dct')
    ret = dctmtx(n).*randsample([1 -1],n,1);
else
    ret = eye(n);
end
