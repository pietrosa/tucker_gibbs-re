function ret = rvmf_gibbs(M, draws)
% Convenience function that performs multiple Gibbs scans to draw from the
%   von Mises-Fisher distribution with a given matrix parameter. An initial
%   draw is made from the uniform distribution on the Stiefel manifold.
%
% Syntax:
%   X = rvmf_gibbs(M, draws) makes draws Gibbs updates to an initial matrix
%       drawn from the uniform distribution on the Stiefel manifold, in order
%       to simulate draws from the von Mises-Fisher distribution with the matrix
%       parameter M.
%
% Inputs:
%   M - [n,p] size, [double] type, matrix parameter for the von Mises Fisher
%       distribution to be simulated
%   draws - [1] size, [int] type, number of updates to be made
%
% Outputs:
%   ret - [n,p] size, [double] type, simulated draw from the von Mises Fisher
%       distribution with matrix parameter M

ret = rustiefel(size(M));
for d = 1:draws
    ret = rvmf_gibbs_update(M, ret);
end
