clc; clear; close all;

path_1 = '/Volumes/Elements/Modularity/Human_10_x2/WSBM_Results';
% streamlines log-normal prior
% path_1 = '/Volumes/Elements/Modularity/Human_10/WSBM_Results';
% % streamlines normal prior

ks = 6:15;

log_evidence = zeros(1,length(ks));
errors = zeros(1,length(ks));

for idx = 1:length(ks)
    k = ks(idx);
    path_2 = sprintf('k_%d', k);
    log_evidence_subj = zeros(1, length(1:10));
    for subj = 1:10
        path_3 = sprintf('subj_%d_k_%d_human_10_LN.mat', subj, k);
        % path_3 = sprintf('subj_%d_k_%d_human_10.mat', subj, k);
        load(fullfile(path_1, path_2, path_3)); 
        % loads in Scores amongst other variables
        log_evidence_subj(subj) = mean(Scores);
    end
    log_evidence(idx) = mean(log_evidence_subj);
    errors(idx) = std(log_evidence_subj);
end

%% Plot

figure;
% plot(ks, log_evidence, 'Linewidth', 2);
e = errorbar(ks, log_evidence, errors/sqrt(10), ...
    '.k', 'MarkerSize', 25, 'LineWidth', 0.1,...
    'MarkerEdgeColor','black','MarkerFaceColor','black');
% e.Color = [0, 0, 0.5430];
% e.LineWidth = 2;
xlabel('Number of Blocks', 'FontSize', 15);
ylabel('Log Likelihood', 'FontSize', 15);
title('Number of Blocks vs. Log Likelihood', 'FontSize', 15);
prettify;