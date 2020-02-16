clc; close all; clear;

W_Distr = 'LogNormal'; % edge weight distribution
E_Distr = 'Bernoulli'; % edge existence distribution
num_iters = 25; % number of WSBM runs

for subj = 1:3
    data_str = sprintf('subj_%d.mat', subj);
    load(data_str);
    A = connectivity; % load -> connectivity, name
    A(A == 0) = NaN; % replace 0s with NaNs
    E = Adj2Edg(A);
    for k = 6:1:15 % sweep over number of communities
        ModelInputs = cell(numel(num_iters),1);
        for iter = 1:num_iters
            ModelInputs{iter} = {k, 'W_Distr', W_Distr, 'E_Distr', E_Distr};
        end
        [Best_Model,Scores,Models] = wsbmLooper(E, ModelInputs);
        filename = sprintf('subj_%d_k_%d_human_10_LN', subj, k);
        save(filename,'Best_Model','Scores','Models');
    end
end

        