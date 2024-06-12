parpool(4);
%% Load hippocampus data
% !!! ADNI data not provided due to privacy restrictions
%load_raw_hippo
% subset variables: radial distance, mTBM1, mTBM2, mTBM3
X = X(:,:,:,1:4,:);

%% Decomposition/compression parameters
R_list = {[5 5 2 4 5], [15 15 2 4 15], [30 30 2 4 30]};
% DR   =  1         0.8               0.6              0.4              0.2
m_list = {size(X), [80 120 2 4 659], [60 90 2 4 494], [40 60 2 4 330], [20 30 2 4 164]};

%% MCMC parameters
ndraws = 1050;
thin = 5;
burnin = 50;
nsave = floor((ndraws-burnin)/thin);
nchains = 4;

%% Run chains and save results
res_all = NaN(nsave,0);
titles_all = [];

%% Main loop
for ir = 1:length(R_list)
    R = R_list{ir}
    for im = 1:length(m_list)
        m = m_list{im}
        setting_err = NaN(nsave,nchains);
        setting_time = NaN(nsave,nchains);
        setting_tau2 = NaN(nsave,nchains);
        draws = {};
        parfor (cc = 1:nchains, 1)
            [draws{cc}, setting_err(:,cc), setting_time(:,cc), ~] = tucker_gibbs(X, R, 'm', m, 'thin', thin, 'draws', ndraws, 'burnin', burnin, 'F_gibbs_draws', [5 5 0 0 5], 'delta2', 1e-1, 'tau2', 1e-1, 'verbose', true);
            setting_tau2(:,cc) = draws{cc}{3};
        end
        pre = strcat("R",string(ir),"_m",string(im));
        titles_all = [titles_all, repelem(strcat(pre,"_err"),nchains), repelem(strcat(pre,"_time"),nchains), repelem(strcat(pre,"_tau2"),nchains)];
        res_all = [res_all, setting_err, setting_time, setting_tau2];
    end
end

writematrix([titles_all;res_all],"hippo_chains.csv");
