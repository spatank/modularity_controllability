clc; close all; clear;

load('data8x3fullscale.mat')

%% Collect edge weights

all_As = A{1,3}; % Scale 125: 234 regions

num_subjects = size(all_As,1);
num_trials = size(all_As,2);

all_raw_data = cell(num_subjects,num_trials);

edge_weights_all = [];

for subj = 1:num_subjects
    for trial = 1:num_trials
        fprintf('Subj = %d, Trial = %d.\n', subj, trial);
        A_raw = all_As{subj,trial}.network; % get adjacency matrix
        A_raw(A_raw == 0) = NaN;
        A_raw(1:size(A_raw,1)+1:end) = NaN; % remove diagonal
        % In the WSBM code and in the QA Networks, missing values are
        % represented by NaNs. However, functions of the BCT typically need
        % 0s for missing values. Hence, conversions of NaNs to 0s are made
        % here.
        A = A_raw;
        edge_list = Adj2Edg(A); 
        edge_weights_all = [edge_weights_all; edge_list(:,3)];
    end
end

%% Plot

figure;
histogram(edge_weights_all);
xlabel('Edge Weight', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);
title('Dataset 1: Skewed Edge Weights', 'FontSize', 15);

figure;
histogram(log(edge_weights_all));
xlabel('Log Edge Weight', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);
title('Dataset 1: Log Edge Weights', 'FontSize', 15);

% figure;
% hist = histfit(log(edge_weights_all));
% hist(1).FaceColor = [.8 .8 1];
% hist(2).Color = [.2 .2 .2];
% xlabel('Log Edge Weight', 'FontSize', 15);
% ylabel('Frequency', 'FontSize', 15);
% title('Dataset 1: Log of Edge Weights', 'FontSize', 15);
% [h, p, stats] = chi2gof(log(edge_weights_all));