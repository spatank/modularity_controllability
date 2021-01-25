clc; close all; clear;

%% Set Network Parameters

dist = 'normal';

N = 234; % number of nodes
density = 0.1485; % desired density
K = round(((N^2-N)*density)/2); % number of edges

Ci = zeros(N,1);
Ci(1:234/2) = 1;
Ci((234/2)+1:end) = 2;

%% Set Controllability Parameters

T = 4; % time horizon for controllability
nor = 0; % matrix normalization flag for stability
thresh = 1; % threshold for slower and faster modes for modal control

noise_flag = 0;

%% Generate Network Ensembles

% iters = 5; % number of networks in each ensemble
% vec = linspace(0.001, 0.95, 5); % number of distinct network topologies

iters = 100; % number of networks in each ensemble
vec = linspace(0.001, 0.95, 25); % number of distinct network topologies

all_nets = zeros(N, N, length(vec), iters);
all_eig_vals = zeros(N, length(vec), iters);
for idx = 1:length(vec)    
    sum_P_1_P_4 = vec(idx);
    P_1 = sum_P_1_P_4/2;
    P_4 = sum_P_1_P_4/2;
    P_2_3 = 1 - sum_P_1_P_4;
    for iter = 1:iters
        A = network_generator(N, K, P_1, P_4, P_2_3, dist, noise_flag);
        all_nets(:,:,idx,iter) = A; % store the adjacency matrix
        all_eig_vals(:,idx,iter) = eig(A); % store eigenvalues
    end
end
max_eig = max(max(max(all_eig_vals))); % max eigenvalue across all ensembles

%% Compute Controllabilities

Q_mat = zeros(length(vec), iters);
avg_ctrb_mat = zeros(length(vec), iters);
mod_ctrb_mat = zeros(length(vec), iters);
min_eng_mat = zeros(length(vec), iters);

for idx = 1:length(vec) 
    fprintf('Idx %d\n', idx);
    for iter = 1:iters
        A = all_nets(:,:,idx,iter);
        Q_mat(idx,iter) = modularity_index(A, Ci);
        A = (A/max_eig);
        avg_ctrb_mat(idx,iter) = mean(avg_ctrb_disc(A, T, nor));
        mod_ctrb_mat(idx,iter) = mean(mod_ctrb_disc(A, 1:size(A,1), thresh, nor));
        min_eng_mat(idx,iter) = mean(min_eng_0_1_node(A, T, 'disc'));
    end
end

Q = mean(Q_mat, 2);
Q_err = std(Q_mat, 0, 2);

avg_ctrb_Z = zscore(avg_ctrb_mat);
avg_ctrb = mean(avg_ctrb_Z, 2);
avg_ctrb_err = std(avg_ctrb_Z,0,2);

mod_ctrb_Z = zscore(mod_ctrb_mat);
mod_ctrb = mean(mod_ctrb_Z, 2);
mod_ctrb_err = std(mod_ctrb_Z,0,2);

min_eng_Z = zscore(min_eng_mat);
min_eng = mean(min_eng_Z, 2);
min_eng_err = std(min_eng_Z,0,2);

%% Plots

path_1 = '/Users/sppatankar/Desktop/Projects/Modularity/Re-submission/';

f = figure('color','w');
errorbar(vec, avg_ctrb, avg_ctrb_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort.: Avg. Ctrb.', 'FontSize', 15);
xlabel('Fraction of Edges in Modules', 'FontSize', 15);
ylabel('Z-Scored Average Controllability', 'FontSize', 15);
prettify
path_2 = strcat('frac_avg_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

figure;
errorbar(vec, mod_ctrb, mod_ctrb_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort: Mod. Ctrb.', 'FontSize', 15);
xlabel('Fraction of Edges in Modules', 'FontSize', 15);
ylabel('Z-Scored Modal Controllability', 'FontSize', 15);
prettify
path_2 = strcat('frac_mod_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

figure;
errorbar(vec, min_eng, min_eng_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort: Min. Eng.', 'FontSize', 15);
xlabel('Fraction of Edges in Modules', 'FontSize', 15);
ylabel('Z-Scored Control Energy', 'FontSize', 15);
prettify
path_2 = strcat('frac_min_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

figure;
errorbar(Q, avg_ctrb, avg_ctrb_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort: Avg. Ctrb.', 'FontSize', 15);
xlabel('Modularity Q', 'FontSize', 15);
ylabel('Z-Scored Average Controllability', 'FontSize', 15);
prettify
path_2 = strcat('Q_avg_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

figure;
errorbar(Q, mod_ctrb, mod_ctrb_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort: Mod. Ctrb.', 'FontSize', 15);
xlabel('Modularity Q', 'FontSize', 15);
ylabel('Z-Scored Modal Controllability', 'FontSize', 15);
prettify
path_2 = strcat('Q_mod_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

figure;
errorbar(Q, min_eng, min_eng_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort.: Min. Eng.', 'FontSize', 15);
xlabel('Modularity Q', 'FontSize', 15);
ylabel('Z-Scored Control Energy', 'FontSize', 15);
prettify
path_2 = strcat('Q_min_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

figure;
errorbar(vec, Q, Q_err, '.k', 'MarkerSize', 25, 'LineWidth', 0.1);
title('Nor.: Disassort. to Assort.: Q and Fraction of Edges in Modules', 'FontSize', 25);
xlabel('Fraction of Edges in Core', 'FontSize', 15);
ylabel('Modularity Q', 'FontSize', 15);
prettify
path_2 = strcat('frac_Q_dis_ass_', dist);
saveas(gcf, fullfile(path_1, path_2), 'epsc')