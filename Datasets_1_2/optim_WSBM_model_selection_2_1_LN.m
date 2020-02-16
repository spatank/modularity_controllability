clc; close all; clear;
addpath(genpath('.'))

%% Run WSBM

load('data8x3fullscale.mat')
all_As = A{1,3}; % Scale 125: 234 regions

% input parameters for WSBM
W_Distr = 'LogNormal';
E_Distr = 'Bernoulli';
num_iters = 10;
num_trials = 3;

for subj = 1:2
    for trial = 1:num_trials
        A = all_As{subj,trial}.network;
        A(1:size(A,1)+1:end) = NaN; % remove diagonal
        A(A == 0) = NaN;
        E = Adj2Edg(A);
        for k = 6:1:15 % sweep over number of communities
            ModelInputs = cell(numel(num_iters),1);
            for iter = 1:num_iters
                ModelInputs{iter} = {k, 'W_Distr', W_Distr, 'E_Distr', E_Distr};
            end
            [Best_Model,Scores,Models] = wsbmLooper(E, ModelInputs);
            filename = sprintf('subj_%d_trial_%d_k_%d_stream_LN', subj, trial, k);
            save(filename,'Best_Model','Scores','Models');
        end
    end
end