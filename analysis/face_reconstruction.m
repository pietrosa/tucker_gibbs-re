%% Load face data
load face_data.mat
X = tensor(X);

%% Random sample of faces
rng(240508)
rind = randsample(400, 16, false);

%% Settings
ndraws = 250;
thin = 10;
burnin = 50;
nsave = floor((ndraws-burnin)/thin);
kk = nsave;
R = [30 30 30];


%% Reconstructions: Effect of compression
tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

% Original data
XX = stack_images(arrayfun(@(k) X.data(:,:,k), rind, 'UniformOutput',false), 4, 4);
nexttile
image(rot90(XX), 'CDataMapping','scaled');
axis off
title("Original")

% R=30: DR=1
[draws, setting_err, setting_time, mixmats] = tucker_gibbs(X, R, 'm', round(size(X).*1), 'thin', thin, 'draws', ndraws, 'burnin', burnin, 'verbose', true);
setting_err
lambda_hat = tensor(squeeze(draws{2}(kk,:,:,:)));
gamma_hat = arrayfun(@(k) squeeze(draws{1}{k}(kk,:,:)), 1:ndims(X), 'UniformOutput', false);
XHAT = ttm(lambda_hat, gamma_hat, 1:ndims(X));
XHAT = ttm(XHAT, mixmats, 1:ndims(X), 't');
XX = stack_images(arrayfun(@(k) XHAT.data(:,:,k), rind, 'UniformOutput',false), 4, 4);
nexttile
image(rot90(XX), 'CDataMapping','scaled');
axis off
title("DR = 1")

% R=30: DR~0.8
[draws, setting_err, setting_time, mixmats] = tucker_gibbs(X, R, 'm', round(size(X).*0.8), 'thin', thin, 'draws', ndraws, 'burnin', burnin, 'verbose', true);
setting_err
lambda_hat = tensor(squeeze(draws{2}(kk,:,:,:)));
gamma_hat = arrayfun(@(k) squeeze(draws{1}{k}(kk,:,:)), 1:ndims(X), 'UniformOutput', false);
XHAT = ttm(lambda_hat, gamma_hat, 1:ndims(X));
XHAT = ttm(XHAT, mixmats, 1:ndims(X), 't');
XX = stack_images(arrayfun(@(k) XHAT.data(:,:,k), rind, 'UniformOutput',false), 4, 4);
nexttile
image(rot90(XX), 'CDataMapping','scaled');
axis off
title("DR = 0.8")

% R=30: DR~0.6
[draws, setting_err, setting_time, mixmats] = tucker_gibbs(X, R, 'm', round(size(X).*0.6), 'thin', thin, 'draws', ndraws, 'burnin', burnin, 'verbose', true);
setting_err
lambda_hat = tensor(squeeze(draws{2}(kk,:,:,:)));
gamma_hat = arrayfun(@(k) squeeze(draws{1}{k}(kk,:,:)), 1:ndims(X), 'UniformOutput', false);
XHAT = ttm(lambda_hat, gamma_hat, 1:ndims(X));
XHAT = ttm(XHAT, mixmats, 1:ndims(X), 't');
XX = stack_images(arrayfun(@(k) XHAT.data(:,:,k), rind, 'UniformOutput',false), 4, 4);
nexttile
image(rot90(XX), 'CDataMapping','scaled');
axis off
title("DR = 0.6")

% R=30: DR~0.4
[draws, setting_err, setting_time, mixmats] = tucker_gibbs(X, R, 'm', round(size(X).*0.8), 'thin', thin, 'draws', ndraws, 'burnin', burnin, 'verbose', true);
setting_err
lambda_hat = tensor(squeeze(draws{2}(kk,:,:,:)));
gamma_hat = arrayfun(@(k) squeeze(draws{1}{k}(kk,:,:)), 1:ndims(X), 'UniformOutput', false);
XHAT = ttm(lambda_hat, gamma_hat, 1:ndims(X));
XHAT = ttm(XHAT, mixmats, 1:ndims(X), 't');
XX = stack_images(arrayfun(@(k) XHAT.data(:,:,k), rind, 'UniformOutput',false), 4, 4);
nexttile
image(rot90(XX), 'CDataMapping','scaled');
axis off
title("DR = 0.4")