function [draws, err, time, mixmats] = tucker_gibbs(X, R, varargin)
% Gibbs sampler for a Bayesian orthogonal Tucker decomposition.
%
% Syntax:
%   [draws, err, time, mixmats] = tucker_gibbs(X, R, varargin) returns saved draws
%       of model components, Frobenius reconstruction error, posterior draw time,
%       and the orthogonal mixing matrices used in an initial preprocessing step.
%
% Inputs:
%   X - [n1,...,nq] size, [tensor] type, input order-q tensor to be decomposed
%   R - [q] size, [int] type, positive ranks for the Tucker decomposition
%   varargin - optional arguments as described below; see main code body for default values
%       - m: [q] size, [int] type, embedding/target dimension for each mode; random
%           embeddings are not applied along modes where m=n
%       - gamma_init: {[n1,R(1)], ..., [nq,R(q)]} size, [double] type, list of
%           initial values for the factor matrices; should be semi-orthogonal
%       - lambda_init: [R(1),...,R(q)] size, [tensor] type, initial value for the
%           core tensor
%       - tau2_init: [1] size, [double] type, positive initial value for the precision
%           parameter tau^2 for X
%       - delta2_init: [1] size, [double] type, positive initial value for the precision
%           parameter delta^2 for the core tensor
%       - F0: {[n1,R(1)],...,[nq,R(q)]} size, [double] type, list of matrix
%           (hyper)parameter for the prior von Mises-Fisher distribution of the
%           factor matrices
%       - lambda0_init: [1] size, [double] type, initial value for the mean of the
%           components of the core tensor
%       - lambda00: [1] size, [double] type, (hyper)parameter for the mean of the
%           prior normal distribution of lambda0
%       - a_delta: [1] size, [double] type, positive (hyper)parameter for the
%           shape of the prior Gamma distribution of delta^2
%       - b_delta: [1] size, [double] type, positive (hyper)parameter for the
%           rate of the prior Gamma distribution of delta^2
%       - a_tau: [1] size, [double] type, positive (hyper)parameter for the
%           shape of the prior Gamma distribution of tau^2
%       - b_tau: [1] size, [double] type, positive (hyper)parameter for the
%           rate of the prior Gamma distribution of tau^2
%       - sigma2: [1] size, [double] type, positive (hyper)parameter for the
%           variance of the prior Normal distribution of lambda0
%       - draws: [1] size, [int] type, total number of posterior updates to be
%           made, including burn-in draws and draws lost due to thinning
%       - burnin: [1] size, [int] type, length of the burn-in period
%       - thin: [1] size, [int] type, thinning factor to apply to the chain; only
%           every thin-th draw is saved
%       - verbose: [1] size, [logical] type, should progress be printed?
%       - use_re_core: [1] size, [logical] type, should random embeddings be used in
%           updates of the core tensor? An experimental option
%       - update_re: [1] size, [logical] type, should the random embeddings be updated
%           every draw (as opposed to being held constant across all draws)? An
%           experimental option
%       - replace: [1] size, [logical] type, should random sampling matrices allow
%           sampling with replacement? An experimental option
%       - mix: [1] size, [logical] type, should the initial mixing step be
%           performed? An experimental option
%       - F_gibbs_draw: [q] size, [int] type, number of Gibbs scans to be be used
%           in the Gibbs subsampler when making posterior draws for the factor
%           matrices; when a component is zero, a rejection sampler (rather than
%           a Gibbs subsampler) is used.
%
% Outputs:
% In the following, nsave = floor((draws-burnin)/thin)
%   draws - list with 6 components containing saved posterior draws for decomposition
%       parameters, as described below
%           - draws{1}: {[nsave,n1,R(1)],...,[nsave,nq,R(q)]} size, [double] type,
%               draws for the factor matrices
%           - draws{2}: [nsave,R(1),...,R(q)] size, [tensor] type, draws for the
%               core tensor
%           - draws{3}: [nsave] size, [double] type, draws for tau2
%           - draws{4}: {[nsave,n1,R(1)],...,[nsave,nq,R(q)]} size, [double] type,
%               matrix parameters for the (von Mises-Fisher) posterior distribution
%               of the factor matrices
%           - draws{5}: [nsave] size, [double] type, draws for lambda0
%           - draws{6}: [nsave] size, [double] type, draws for delta2
%   err - [nsave] size, [double] type, Frobenius reconstruction error computed
%       at every draw where posterior draws were saved
%   time - [nsave] size, [double] type, total time for all posterior draws, computed
%       for every draw where posterior draws were saved
%   mixmats - {[n1,n1],...,[nq,nq]}, [double] type, orthogonal mixing matrices
%       applied in an initial preprocessing step. The transpose of these should
%       be multipled into draws for the factor matrices in order to revert back
%       to the original tensor space

%% extract properties of input tensor
q = ndims(X);
n = size(X);
modes = 1:q;
rust = arrayfun(@(k) rustiefel([n(k),R(k)]), 1:q, 'UniformOutput', false);


%% algorithm parameters
params = inputParser;
% no check for gamma_init, lambda_init, or F0 currently implemented
params.addParameter('m', n, @(x) length(x)==q && all(x<=n));
params.addParameter('gamma_init', rust, @(x) 1);
params.addParameter('lambda_init', ttm(X, rust, 1:q, 't'), @(x) 1);
params.addParameter('tau2_init', 10e-4, @(x) isscalar(x) && x>0);
params.addParameter('delta2_init', 10e-4, @(x) isscalar(x) && x>0);
params.addParameter('F0', arrayfun(@(j) zeros(n(j),R(j)), 1:q, 'UniformOutput', false), @(x) 1);
params.addParameter('lambda0_init', 0, @isscalar);
params.addParameter('lambda00', 0, @isscalar);
params.addParameter('a_delta', 10e-4, @(x) isscalar(x) && x>0);
params.addParameter('b_delta', 10e-4, @(x) isscalar(x) && x>0);
params.addParameter('a_tau', 10e-4, @(x) isscalar(x) && x>0);
params.addParameter('b_tau', 10e-4, @(x) isscalar(x) && x>0);
params.addParameter('sigma2', 10, @(x) isscalar(x) && x>0);
params.addParameter('draws', 10, @(x) isscalar(x) && x>=0 && floor(x)==x);
params.addParameter('burnin', 0, @(x) isscalar(x) && x>=0 && floor(x)==x);
params.addParameter('thin', 10, @(x) isscalar(x) && x>=0 && floor(x)==x);
params.addParameter('verbose', true, @(x) isscalar(x) && islogical(x));
params.addParameter('use_re_core', true, @(x) isscalar(x) && islogical(x));
params.addParameter('update_re', true, @(x) isscalar(x) && islogical(x));
params.addParameter('replace', false, @(x) isscalar(x) && islogical(x));
params.addParameter('mix', true, @(x) isscalar(x) && islogical(x));
params.addParameter('F_gibbs_draws', repelem(5,q), @(x) length(x)==q);

params.parse(varargin{:});

%% copy parameters
gamma = params.Results.gamma_init;
lambda = params.Results.lambda_init;
tau2 = params.Results.tau2_init;
delta2 = params.Results.delta2_init;
F0 = params.Results.F0;
lambda0 = params.Results.lambda0_init;
lambda00 = params.Results.lambda00;
a_delta = params.Results.a_delta;
b_delta = params.Results.b_delta;
a_tau = params.Results.a_tau;
b_tau = params.Results.b_tau;
sigma2 = params.Results.sigma2;
draws = params.Results.draws;
burnin = params.Results.burnin;
thin = params.Results.thin;
verbose = params.Results.verbose;
use_re_core = params.Results.use_re_core;
update_re = params.Results.update_re;
replace = params.Results.replace;
mix = params.Results.mix;
m = params.Results.m;
F_gibbs_draws = params.Results.F_gibbs_draws;

%% apply mixing
if mix
    mixmats = arrayfun(@(k) mixmat(n(k), 'dct'), modes, 'UniformOutput', false);
    X = ttm(X, mixmats, modes);
else
    mixmats = NaN;
end

%% calculate other quantities
nn = prod(n);
RR = prod(R);
mm = prod(m);
use_re = any(m~=n);
if use_re
    a = m./n;
    b = arrayfun(@(j) prod(a(setdiff(modes,j))), modes);
    b0 = prod(a);
end

%% save draws
nsave = floor((draws-burnin)/thin);
gamma_draws = arrayfun(@(k) NaN(nsave,n(k),R(k)), modes, 'UniformOutput', false);
lambda_draws = NaN([nsave,R]);
Rcolon = repmat({':'},1,q);
tau2_draws = NaN(nsave,1);
F_draws = arrayfun(@(k) NaN(nsave,n(k),R(k)), modes, 'UniformOutput', false);
lambda0_draws = NaN(nsave,1);
delta2_draws = NaN(nsave,1);
err = NaN(nsave,1);
time = NaN(nsave,1);

%% some preallocation
vecl = arrayfun(@(k) NaN(R(k)), modes, 'UniformOutput', false);
diagl = arrayfun(@(k) NaN(R(k),1), modes, 'UniformOutput', false);
F_post = arrayfun(@(k) NaN(n(k),R(k)), modes, 'UniformOutput', false);

%% main loop
for s = 1:draws
    T0 = tic;
    % create subsample indices
    if use_re && (update_re || s==1)
        sinds = arrayfun(@(k) randsample(n(k),m(k),replace), modes, 'UniformOutput', false);
    end

    % update factor matrices
    for j = modes
        modes_ = setdiff(modes, j);
        if use_re
            sinds_ = listreplace(sinds, j, ':');
            X_ = X(sinds_{:});
            gamma_ = arrayfun(@(k) gamma{k}(sinds{k},:), modes_, 'UniformOutput', false);
            % obtain posterior
            tmp = ttm(X_, gamma_, modes_, 't');
            F_post{j} = b(j)*tau2 * tenmat(tmp, j).data * tenmat(lambda, j).data' + F0{j};
        else
            tmp = ttm(X, gamma(modes_), modes_, 't');
            F_post{j} = tau2 * tenmat(tmp, j).data * tenmat(lambda, j).data' + F0{j};
        end

        % draw from posterior
        if F_gibbs_draws(j)
            gamma{j} = rvmf_gibbs(F_post{j}, F_gibbs_draws(j));
        else
            gamma{j} = rvmf_reject(F_post{j});
        end
    end

    % update core tensor
    if use_re && use_re_core
        X_ = X(sinds{:});
        gamma_ = arrayfun(@(k) gamma{k}(sinds{k},:), modes, 'UniformOutput', false);
        for j = modes
            [V,D] = eig(xtx(gamma_{j}));
            vecl{j} = V;
            diagl{j} = diag(D)';
        end
        diagl_flip = flip(diagl);
        Dinv = 1./(b0*tau2 * kron_diag(diagl_flip{:}) + delta2)';
        tAGV = arrayfun(@(k) (gamma_{k} * vecl{k})', modes, 'UniformOutput', false);
        V_colsum_rev = arrayfun(@(k) sum(vecl{k},1), flip(modes), 'UniformOutput', false);
        tmp = ttm(X_, tAGV, modes);
        lambda_post_mean = Dinv .* (b0*tau2 * vec(tmp.data) + lambda0*delta2 * kron_diag(V_colsum_rev{:})');
        lambda(:) = lambda_post_mean + randn(RR,1).*sqrt(Dinv);
        lambda = ttm(lambda, vecl, modes);
    else
        lambda_post_var = 1/(tau2+delta2);
        lambda_post_mean = lambda_post_var * (tau2 * vec(ttm(X, gamma, modes, 't').data) + delta2*lambda0);
        lambda(:) = lambda_post_mean + randn(RR, 1)*sqrt(lambda_post_var);
    end

    % update precisions
    if use_re && use_re_core
        shape = a_tau + 0.5*mm;
        rate =  b_tau + 0.5*b0*norm(X_ - ttm(lambda, gamma_, modes))^2;
    else
        shape = a_tau + 0.5*nn;
        rate = b_tau + 0.5*norm(X - ttm(lambda, gamma, modes))^2;
    end
    tau2 = gamrnd(shape, 1/rate);

    shape = a_delta + 0.5*RR;
    rate = b_delta + 0.5*norm(lambda - lambda0)^2;
    delta2 = gamrnd(shape, 1/rate);

    % update core mean
    lambda0_var = sigma2 / (1+ sigma2*delta2*RR);
    lambda0_mean = lambda0_var * (lambda00/sigma2 + delta2*sum(lambda.data, 'all'));
    lambda0 = lambda0_mean + randn()*sqrt(lambda0_var);

    % save draws
    if s>burnin && mod(s-burnin, thin)==0
        ss = (s-burnin)/thin;
        % time for update
        time(ss) = toc(T0);
        % draws
        for j = modes
            gamma_draws{j}(ss,:,:) = gamma{j};
            F_draws{j}(ss,:,:) = F_post{j};

        end
        ssRcolon = [{ss},Rcolon(:)'];
        lambda_draws(ssRcolon{:}) = lambda.data;
        if use_re
            % to account for not every applying the scale factors b(j) or b0 to X_
            tau2_draws(ss) = tau2*b0;
        else
            tau2_draws(ss) = tau2;
        end
        lambda0_draws(ss) = lambda0;
        delta2_draws(ss) = delta2;
        % error
        Xhat = ttm(lambda, gamma, modes);
        err(ss) = norm(X-Xhat);
    end

    if verbose
        if s<=burnin
            fprintf('(burnin) ');
        end
        if s>burnin && mod(s-burnin, thin)==0
            fprintf('iter = %i\ttime = %.2f\terr = %.4f\n', s, time(ss), err(ss));
        else
            fprintf('iter = %i\ttime = %.2f\n', s, toc(T0));
        end
    end
end
draws = {gamma_draws, lambda_draws, tau2_draws, F_draws, lambda0_draws, delta2_draws};
end
