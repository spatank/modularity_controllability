clc; clear; close all;

% path_1 = ['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
%     'Human_8_3/Data_Scripts/WSBM_Results']; % QA with normal prior
% path_1 = ['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
%     'Human_8_3/Data_Scripts_Streamlines_LogNormal/WSBM_Results_Streamlines_LN']; 
% streamlines with log-normal prior
path_1 = ['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
    'Human_8_3/Data_Scripts_Streamlines_Normal/WSBM_Results_Streamlines_Normal']; 
% streamlines with normal prior

ks = 6:15;
num_subjects = 8;
num_trials = 3;

log_evidence = zeros(num_trials, length(ks));
errors = zeros(num_trials, length(ks));

for trial = 1:num_trials
    cd(path_1); % restore path to root of WSBM_Results
    path_2 = sprintf('Trial %d', trial);
    for idx = 1:length(ks)
        k = ks(idx);
        path_3 = sprintf('k_%d', k);
        log_evidence_subj = zeros(1, num_subjects);
        for subj = 1:num_subjects
%             path_4 = sprintf('subj_%d_trial_%d_k_%d.mat', ...
%                 subj, trial, k); % QA Gaussian
%             path_4 = sprintf('subj_%d_trial_%d_k_%d_stream_LN.mat', ...
%                 subj, trial, k); % stream LogNormal
            path_4 = sprintf('subj_%d_trial_%d_k_%d_stream.mat', ...
                subj, trial, k); % stream Normal
            % load all scores for current trial at the fixed k
            load(fullfile(path_1, path_2, path_3, path_4)); 
            % loads in Scores amongst other variables
            log_evidence_subj(subj) = mean(Scores);
        end
        log_evidence(trial, idx) = mean(log_evidence_subj);
        errors(trial, idx) = std(log_evidence_subj);
    end
end


%% Plot

figure;
e = errorbar(ks, mean(log_evidence), std(log_evidence, 0, 1)/sqrt(8),...
    '.k', 'MarkerSize', 25, 'LineWidth', 0.1, ...
    'MarkerEdgeColor','black','MarkerFaceColor','black');
xlabel('Number of Blocks', 'FontSize', 15);
ylabel('Log Likelihood', 'FontSize', 15);
title('Number of Blocks vs. Log Likelihood', 'FontSize', 15);
prettify