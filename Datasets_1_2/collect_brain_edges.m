clc; close all; clear;

load(['/Users/sppatankar/Developer/modularity_controllability/Data/', ...
    'Human_8_3/ncomm_8_3_qa_gfa.mat']);

%% Collect edge weights

% all_As = qaNetworks(:,:,1); % Scale 33: 83 regions
% all_As = qaNetworks(:,:,2); % Scale 60: 129 regions
all_As = qaNetworks(:,:,3); % Scale 125: 234 regions
% all_As = qaNetworks(:,:,4); % Scale 250: 463 regions
% all_As = qaNetworks(:,:,5); % Scale 500: 1015 regions

num_subjects = size(all_As,1);
num_trials = size(all_As,2);

all_raw_data = cell(num_subjects,num_trials);

edge_weights_all = [];

for subj = 1:num_subjects
    for trial = 1:num_trials
        fprintf('Subj = %d, Trial = %d.\n', subj, trial);
        A_raw = all_As{subj,trial}; % get adjacency matrix
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
hist = histfit(edge_weights_all);
hist(1).FaceColor = [.8 .8 1];
hist(2).Color = [.2 .2 .2];
xlabel('Edge Weight', 'FontSize', 15);
ylabel('Frequency', 'FontSize', 15);
title('Dataset 1: Normally Distributed Edge Weights', 'FontSize', 15);

[h, p, stats] = chi2gof(edge_weights_all);