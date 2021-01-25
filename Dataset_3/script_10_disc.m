clc; close all; clear;

data_set = 'd_3_1';

load('/Volumes/My Passport/Modularity_2/Human_8_3/yeo_atlas.mat')
c = yeo_atlas;
% for reference:
% VIS = 1;
% SOM = 2;
% DOR = 3;
% VEN = 4;
% LIM = 5;
% FPC = 6;
% DMN = 7;
% SUB = 8;

%% Generate Data Structure

% input parameters for controllability metrics
T_avg = 4; % time horizon for average controllability and minimum control energy
T_eng = 4; % time horizon for minimum control energy
nor = 0; % matrix normalization flag for stability
thresh = 1; % threshold for slower and faster modes for modal control

if strcmp(data_set, 'd_3_1')
    path_1 = '/Volumes/My Passport/Modularity_2/Human_10_x2/WSBM_Results/'; % LN
    k = 10; % number of communities; with LogNormal prior
end
if strcmp(data_set, 'd_3_2')
    path_1 = '/Volumes/My Passport/Modularity_2/Human_10/WSBM_Results/'; % Normal
    k = 11; % number of communities; with Normal prior
end

num_subjects = 10;

all_raw_data = cell(1, num_subjects);

for subj = 1:num_subjects

    subj_struct = struct('A', [], 'M', [], ...
    'part_coeff', [], 'within_mod_z', [], ...
    'node_strength', [], 'sc', [], ...
    'avg_ctrb_disc',[], 'mod_ctrb_disc',[], ...
    'min_eng', []);

    fprintf('Subj = %d.\n', subj);
    path_2 = sprintf('Data_Copy/Subj_%d.mat', subj);
    load(fullfile(path_1, path_2)); % loads 'connectivity' and 'name'
    
    A_raw = connectivity;
    
    A = (A_raw/(eigs(A_raw,1)));
    
    subj_struct.A = A;
    
    if strcmp(data_set, 'd_3_1')
        path_3 = sprintf('k_%d/subj_%d_k_%d_human_10_LN.mat', k, subj, k);
    end
    if strcmp(data_set, 'd_3_2')
        path_3 = sprintf('k_%d/subj_%d_k_%d_human_10.mat', k, subj, k);
    end

    load(fullfile(path_1, path_3)); % loads 'Best_Model', 'Models', and 'Scores'

    Ci = zeros(size(A, 1), length(Models)); % multiple fitted WSBMs 
    for idx = 1:length(Models)
        model = Models(idx);
        labels = assign_communities(model{1, 1});
        Ci(:, idx) = labels;
    end
    
    M = central_partition(Ci);
    subj_struct.M = M;

    subj_struct.part_coeff = participation_coef(A, M);
    subj_struct.within_mod_z = module_degree_zscore(A, M);
    subj_struct.node_strength = strengths_und(A);

    subj_struct.spec_met = diag(expm(A));
    
    subj_struct.avg_ctrb_disc = avg_ctrb_disc(A, T_avg, nor);
    subj_struct.mod_ctrb_disc = mod_ctrb_disc(A, 1:size(A, 1), thresh, nor); 
    
    subj_struct.min_eng = min_eng_0_1_node(A, T_eng, 'disc'); 
    
    all_raw_data{subj} = subj_struct;
end

%% Average the Values Across Subjects

part_coeff_all = [];
within_mod_z_all = [];
node_strength_all = [];
sc_all = [];
avg_ctrb_all = [];
mod_ctrb_all = [];
min_eng_all = [];

for subj = 1:num_subjects
    part_coeff_all(:, subj) = all_raw_data{1, subj}.part_coeff;
    within_mod_z_all(:, subj) = all_raw_data{1, subj}.within_mod_z;
    node_strength_all(:, subj) = all_raw_data{1, subj}.node_strength;
    sc_all(:, subj) = all_raw_data{1, subj}.spec_met;
    avg_ctrb_all(:, subj) = all_raw_data{1, subj}.avg_ctrb_disc;
    mod_ctrb_all(:, subj) = all_raw_data{1, subj}.mod_ctrb_disc;
    min_eng_all(:, subj) = all_raw_data{1, subj}.min_eng;
end

part_coeff_mean = mean(part_coeff_all, 2);
within_mod_z_mean = mean(within_mod_z_all, 2);
node_strength_mean = mean(node_strength_all, 2);
sc_mean = mean(sc_all, 2);

avg_ctrb_mean = mean(avg_ctrb_all, 2);
mod_ctrb_mean = mean(mod_ctrb_all, 2);
min_eng_mean = mean(min_eng_all, 2);


%% Correlations with Degree 

path_1 = '/Users/sppatankar/Desktop/Projects/Modularity/Re-submission/';

fileID = fopen(fullfile(path_1, sprintf('%s.txt', data_set)),'w');

[r_NS_PC, p_NS_PC] = corr(node_strength_mean, part_coeff_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. part_coeff ~ node_strength\n', r_NS_PC, p_NS_PC);
[r_NS_Z, p_NS_Z] = corr(node_strength_mean, within_mod_z_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. within_mod_z ~ node_strength\n', r_NS_Z, p_NS_Z);
[r_NS_SC, p_NS_SC] = corr(node_strength_mean, sc_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. sc ~ node_strength\n', r_NS_SC, p_NS_SC);

fprintf(fileID, '\n');

%% Minimum Control Energy Models

% Univariate Models

[r_NS_MIN, p_NS_MIN] = corr(node_strength_mean, min_eng_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ node_strength\n', r_NS_MIN, p_NS_MIN);
[r_PC_MIN, p_PC_MIN] = corr(part_coeff_mean, min_eng_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ part_coeff\n', r_PC_MIN, p_PC_MIN);
[r_Z_MIN, p_Z_MIN] = corr(within_mod_z_mean, min_eng_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ within_mod_deg_z\n', r_Z_MIN, p_Z_MIN);
[r_SC_MIN, p_SC_MIN] = corr(sc_mean, min_eng_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ sc\n', r_SC_MIN, p_SC_MIN);

% Multivariate Models

[r_PC_NS_MIN, p_PC_NS_MIN, residuals] = ... 
    partialcorr_with_resids(part_coeff_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ part_coeff, controlling for node_strength\n', r_PC_NS_MIN, p_PC_NS_MIN);
PC_resid = residuals(:, 1);
MIN_resid = residuals(:, 2);
coeffs = polyfit(PC_resid, MIN_resid, 1);
x = linspace(min(PC_resid), max(PC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(PC_resid, MIN_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Participation Coefficient Residual', 'FontSize', 15);
ylabel('Minimum Control Energy Residual', 'FontSize', 15);
title('Participation Coefficient vs. Minimum Control Energy', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('PC_MIN_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

[r_Z_NS_MIN, p_Z_NS_MIN, residuals] = ... 
    partialcorr_with_resids(within_mod_z_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ within_mod_deg_z, controlling for node_strength\n', r_Z_NS_MIN, p_Z_NS_MIN);
Z_resid = residuals(:, 1);
MIN_resid = residuals(:, 2);
coeffs = polyfit(Z_resid, MIN_resid, 1);
x = linspace(min(Z_resid), max(Z_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(Z_resid, MIN_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Within-Module Degree Z-Score Residual', 'FontSize', 15);
ylabel('Minimum Control Energy Residual', 'FontSize', 15);
title('Within-Module Degree Z-Score vs. Minimum Control Energy', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('Z_MIN_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

[r_SC_NS_MIN, p_SC_NS_MIN, residuals] = ... 
    partialcorr_with_resids(sc_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ sc, controlling for node_strength\n', r_SC_NS_MIN, p_SC_NS_MIN);
SC_resid = residuals(:, 1);
MIN_resid = residuals(:, 2);
coeffs = polyfit(SC_resid, MIN_resid, 1);
x = linspace(min(SC_resid), max(SC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SC_resid, MIN_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Subgraph Centrality Residual','FontSize',15);
ylabel('Minimum Control Energy Residual','FontSize',15);
title('SC vs. Minimum Control Energy', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('SC_MIN_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

fprintf(fileID, '\n');

%% Average Controllability Models

% Univariate Models

[r_NS_AVG, p_NS_AVG] = corr(node_strength_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ node_strength\n', r_NS_AVG, p_NS_AVG);
[r_PC_AVG, p_PC_AVG] = corr(part_coeff_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ part_coeff\n', r_PC_AVG, p_PC_AVG);
[r_Z_AVG, p_Z_AVG] = corr(within_mod_z_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ within_mod_deg_z\n', r_Z_AVG, p_Z_AVG);
[r_SC_AVG, p_SC_AVG] = corr(sc_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ sc\n', r_SC_AVG, p_SC_AVG);

% Multivariate Models

[r_PC_NS_AVG, p_PC_NS_AVG, residuals] = ... 
    partialcorr_with_resids(part_coeff_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ part_coeff, controlling for node_strength\n', r_PC_NS_AVG, p_PC_NS_AVG);
PC_resid = residuals(:, 1);
AVG_resid = residuals(:, 2);
coeffs = polyfit(PC_resid, AVG_resid, 1);
x = linspace(min(PC_resid), max(PC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(PC_resid, AVG_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Participation Coefficient Residual', 'FontSize', 15);
ylabel('Average Controllability Residual', 'FontSize', 15);
title('Participation Coefficient vs. Average Controllability',...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('PC_AVG_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

[r_Z_NS_AVG, p_Z_NS_AVG, residuals] = ... 
    partialcorr_with_resids(within_mod_z_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ within_mod_deg_z, controlling for node_strength\n', r_Z_NS_AVG, p_Z_NS_AVG);
Z_resid = residuals(:, 1);
AVG_resid = residuals(:, 2);
coeffs = polyfit(Z_resid, AVG_resid, 1);
x = linspace(min(Z_resid), max(Z_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(Z_resid, AVG_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Within-Module Degree Z-Score Residual', 'FontSize', 15);
ylabel('Average Controllability Residual', 'FontSize', 15);
title('Within-Module Degree Z-Score vs. Average Controllability', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('Z_AVG_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

[r_SC_NS_AVG, p_SC_NS_AVG, residuals] = ... 
    partialcorr_with_resids(sc_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ sc, controlling for node_strength\n', r_SC_NS_AVG, p_SC_NS_AVG);
SC_resid = residuals(:, 1);
AVG_resid = residuals(:, 2);
coeffs = polyfit(SC_resid, AVG_resid, 1);
x = linspace(min(SC_resid), max(SC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SC_resid, AVG_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Subgraph Centrality Residual', 'FontSize', 15);
ylabel('Average Controllability Residual', 'FontSize', 15);
title('SC vs. Average Controllability', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('SC_AVG_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')


fprintf(fileID, '\n');

%% Modal Controllability Models

% Univariate Models

[r_NS_MOD, p_NS_MOD] = corr(node_strength_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ node_strength\n', r_NS_MOD, p_NS_MOD);
[r_PC_MOD, p_PC_MOD] = corr(part_coeff_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ part_coeff\n', r_PC_MOD, p_PC_MOD);
[r_Z_MOD, p_Z_MOD] = corr(within_mod_z_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ within_mod_deg_z\n', r_Z_MOD, p_Z_MOD);
[r_SC_MOD, p_SC_MOD] = corr(sc_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ sc\n', r_SC_MOD, p_SC_MOD);

% Multivariate Models

[r_PC_NS_MOD, p_PC_NS_MOD, residuals] = ... 
    partialcorr_with_resids(part_coeff_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ part_coeff, controlling for node_strength\n', r_PC_NS_MOD, p_PC_NS_MOD);
PC_resid = residuals(:, 1);
MOD_resid = residuals(:, 2);
coeffs = polyfit(PC_resid, MOD_resid, 1);
x = linspace(min(PC_resid), max(PC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(PC_resid, MOD_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Participation Coefficient Residual', 'FontSize', 15);
ylabel('Modal Controllability Residual', 'FontSize', 15);
title('Participation Coefficient vs. Modal Controllability', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('PC_MOD_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

[r_Z_NS_MOD, p_Z_NS_MOD, residuals] = ... 
    partialcorr_with_resids(within_mod_z_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ within_mod_deg_z, controlling for node_strength\n', r_Z_NS_MOD, p_Z_NS_MOD);
Z_resid = residuals(:, 1);
MOD_resid = residuals(:, 2);
coeffs = polyfit(Z_resid, MOD_resid, 1);
x = linspace(min(Z_resid), max(Z_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(Z_resid, MOD_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Within-Module Degree Z-Score Residual', 'FontSize', 15);
ylabel('Modal Controllability Residual', 'FontSize', 15);
title('Within-Module Degree Z-Score vs. Modal Controllability', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('Z_MOD_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

[r_SC_NS_MOD, p_SPEC_NS_MOD, residuals] = ... 
    partialcorr_with_resids(sc_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf(fileID, 'r = %f, p = %f. mod_ctrb ~ sc, controlling for node_strength\n', r_SC_NS_MOD, p_SPEC_NS_MOD);
SC_resid = residuals(:, 1);
MOD_resid = residuals(:, 2);
coeffs = polyfit(SC_resid, MOD_resid, 1);
x = linspace(min(SC_resid), max(SC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SC_resid, MOD_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Subgraph Centrality Residual', 'FontSize', 15);
ylabel('Modal Controllability Residual', 'FontSize', 15);
title('SC vs. Modal Controllability', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('SC_MOD_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')