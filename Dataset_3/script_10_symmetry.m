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
    'node_strength', [], 'sc', [], ...
    'avg_ctrb_disc',[], 'min_eng', []);

    fprintf('Subj = %d.\n', subj);
    path_2 = sprintf('Data_Copy/Subj_%d.mat', subj);
    load(fullfile(path_1, path_2)); % loads 'connectivity' and 'name'
    
    A_raw = connectivity;
    
    A_raw = reduce_edge_strength(A_raw, 50, 50); % reduce 

    A = (A_raw/(eigs(A_raw,1)));
    
    subj_struct.A = A;
    
    if strcmp(data_set, 'd_3_1')
        path_3 = sprintf('k_%d/subj_%d_k_%d_human_10_LN.mat', k, subj, k);
    end
    if strcmp(data_set, 'd_3_2')
        path_3 = sprintf('k_%d/subj_%d_k_%d_human_10.mat', k, subj, k);
    end

    load(fullfile(path_1, path_3)); % loads 'Best_Model', 'Models', and 'Scores'

    subj_struct.node_strength = strengths_und(A);

    subj_struct.spec_met = diag(expm(A));
    
    subj_struct.avg_ctrb_disc = avg_ctrb_disc(A, T_avg, nor);
    
    subj_struct.min_eng = min_eng_0_1_node(A, T_eng, 'disc'); 
    
    all_raw_data{subj} = subj_struct;
end

%% Average the Values Across Subjects

node_strength_all = [];
sc_all = [];
avg_ctrb_all = [];
min_eng_all = [];

for subj = 1:num_subjects
    node_strength_all(:, subj) = all_raw_data{1, subj}.node_strength;
    sc_all(:, subj) = all_raw_data{1, subj}.spec_met;
    avg_ctrb_all(:, subj) = all_raw_data{1, subj}.avg_ctrb_disc;
    min_eng_all(:, subj) = all_raw_data{1, subj}.min_eng;
end

node_strength_mean = mean(node_strength_all, 2);
sc_mean = mean(sc_all, 2);

avg_ctrb_mean = mean(avg_ctrb_all, 2);
min_eng_mean = mean(min_eng_all, 2);


%% Correlations with Degree 

path_1 = '/Users/sppatankar/Desktop/Projects/Modularity/Re-submission/symmetry/';

fileID = fopen(fullfile(path_1, sprintf('%s_asym.txt', data_set)),'w');

[r_NS_SC, p_NS_SC] = corr(node_strength_mean, sc_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. sc ~ node_strength\n', r_NS_SC, p_NS_SC);

fprintf(fileID, '\n');

%% Minimum Control Energy Models

% Univariate Models

[r_NS_MIN, p_NS_MIN] = corr(node_strength_mean, min_eng_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ node_strength\n', r_NS_MIN, p_NS_MIN);
[r_SC_MIN, p_SC_MIN] = corr(sc_mean, min_eng_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. min_eng ~ sc\n', r_SC_MIN, p_SC_MIN);

% Multivariate Model

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
path_2 = strcat('asym_SC_MIN_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

fprintf(fileID, '\n');

%% Average Controllability Models

% Univariate Models

[r_NS_AVG, p_NS_AVG] = corr(node_strength_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ node_strength\n', r_NS_AVG, p_NS_AVG);
[r_SC_AVG, p_SC_AVG] = corr(sc_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf(fileID, 'r = %f, p = %f. avg_ctrb ~ sc\n', r_SC_AVG, p_SC_AVG);

% Multivariate Model

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
path_2 = strcat('asym_SC_AVG_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')


fprintf(fileID, '\n');