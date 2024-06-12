parpool(4);
%% Load face data
load face_data.mat
X = tensor(X);

%% Decomposition/compression parameters
R_list = {[5 5 5], [15 15 15], [30 30 30]};
% DR   =  1         0.8          0.6          0.4          0.2
m_list = {size(X), [74 90 320], [55 67 240], [37 45 160], [18 22 80]};

%% MCMC parameters
ndraws = 550;
thin = 5;
burnin = 50;
nsave = floor((ndraws-burnin)/thin);
nchains = 4;

%% Run chains and save results
res_all = NaN(nsave,0);
titles_all = [];

for ir = 1:length(R_list)
    R = R_list{ir}
    for im = 1:length(m_list)
        m = m_list{im}
        setting_err = NaN(nsave, nchains);
        setting_time = NaN(nsave, nchains);
        setting_tau2 = NaN(nsave, nchains);
        draws = {};
        parfor (cc = 1:nchains,4)
            [draws{cc}, setting_err(:,cc), setting_time(:,cc), ~] = tucker_gibbs(X, R, 'm', m, 'thin', thin, 'draws', ndraws, 'burnin', burnin, 'verbose', true);
            setting_tau2(:,cc) = draws{cc}{3};
        end
        pre = strcat("R",string(ir),"_m",string(im));
        titles_all = [titles_all, repelem(strcat(pre,"_err"),nchains), repelem(strcat(pre,"_time"),nchains), repelem(strcat(pre,"_tau2"),nchains)];
        res_all = [res_all, setting_err, setting_time, setting_tau2];
    end
end

writematrix([titles_all;res_all],"face_chains.csv");
