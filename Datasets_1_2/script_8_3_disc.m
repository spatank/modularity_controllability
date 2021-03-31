clc; clear; close all;

addpath(genpath(['/Users/sppatankar/Developer/modularity_controllability/'...
    'Helper']))

data_set = 'd_2_1';

if strcmp(data_set, 'd_1')
    load(['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
        'Human_8_3/ncomm_8_3_qa_gfa.mat']); % for QA weighted streamlines
end
if strcmp(data_set, 'd_2_1') || strcmp(data_set, 'd_2_2')
    load(['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
        'Human_8_3/data8x3fullscale.mat']); % for streamlines
end

%% Generate Data Structure

% input parameters for controllability metrics
T_avg = 4; % time horizon for average controllability 
T_eng = 4; % time horizon for minimum control energy
thresh = 1; % threshold for slower and faster modes for modal control

if strcmp(data_set, 'd_1') % QA with normal prior
    all_As = qaNetworks(:,:,3); % Scale 125: 234 regions
    path_1 = ['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
        'Human_8_3/Data_Scripts/WSBM_Results']; 
    k = 12; % number of communities for QA  with Gaussian prior
end
if strcmp(data_set, 'd_2_1') % streamlines with log-normal prior
    all_As = A{1,3}; % Scale 125: 234 regions
    path_1 = ['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
        'Human_8_3/Data_Scripts_Streamlines_LogNormal/WSBM_Results_Streamlines_LN'];
    k = 14; % number of communities for streamlines with LN prior
end
if strcmp(data_set, 'd_2_2') % streamlines with Gaussian prior
    all_As = A{1,3}; % Scale 125: 234 regions
    path_1 = ['/Users/sppatankar/Developer/modularity_controllability/Data/' ...
        'Human_8_3/Data_Scripts_Streamlines_Normal/WSBM_Results_Streamlines_Normal']; 
    k = 12; % number of communities for streamlines with Gaussian prior
end

num_subjects = size(all_As, 1);
num_trials = size(all_As, 2);

all_raw_data = cell(num_subjects,num_trials);

for subj = 1:num_subjects

    subj_struct = struct('A',[],'M',[],...
    'part_coeff',[],'within_mod_z',[],...
    'node_strength',[],'sc',[],...
    'avg_ctrb_disc',[],'mod_ctrb_disc',[],...
    'min_eng',[]);
    
    for trial = 1:num_trials
        
        path_2 = sprintf('Trial %d', trial);
        
        fprintf('Subj = %d, Trial = %d.\n', subj, trial);
        
        if strcmp(data_set, 'd_1')
            A_raw = all_As{subj,trial}; % get adjacency matrix for QA weighted networks
            path_3 = sprintf('k_%d/subj_%d_trial_%d_k_%d', k, subj, trial, k);
        end
        if strcmp(data_set, 'd_2_1')
            A_raw = all_As{subj,trial}.network; % get adjacency matrix for streamlines
            path_3 = sprintf('k_%d/subj_%d_trial_%d_k_%d_stream_LN', k, subj, trial, k);
        end
        if strcmp(data_set, 'd_2_2')
            A_raw = all_As{subj,trial}.network; % get adjacency matrix for streamlines
            path_3 = sprintf('k_%d/subj_%d_trial_%d_k_%d_stream', k, subj, trial, k);
        end

        A_raw(1:size(A_raw,1)+1:end) = NaN; % remove diagonal
        % In the WSBM code and in the QA Networks, missing values are
        % represented by NaNs. However, functions of the BCT typically need
        % 0s for missing values. Hence, conversions of NaNs to 0s are made
        % here.
        A_raw(isnan(A_raw)) = 0;
             
        A = (A_raw/(eigs(A_raw, 1)));
        
        subj_struct.A = A;

        load(fullfile(path_1, path_2, path_3)); % loads 'Best_Model', 'Models', and 'Scores'
        
        M = central_partition(Models);
        subj_struct.M = M;
        
        subj_struct.part_coeff = participation_coef(A, M);
        subj_struct.within_mod_z = module_degree_zscore(A, M);
        subj_struct.node_strength = strengths_und(A);
 
        subj_struct.sc = diag(expm(A)); % communicability
        
        subj_struct.avg_ctrb_disc = avg_ctrb_disc(A, T_avg);
        subj_struct.mod_ctrb_disc = mod_ctrb_disc(A, 1:size(A,1), thresh);

        subj_struct.min_eng = min_eng_0_1_node(A, T_eng, 'disc');

        all_raw_data{subj,trial} = subj_struct;
    end
end


%% Average the Values Across the Three Scanning Sessions

data_struct = struct('A', [], ...
    'part_coeff', [], 'within_mod_z', [], ...
    'node_strength', [],'sc', [], ...
    'avg_ctrb_disc', [], 'mod_ctrb_disc', [], ...
    'min_eng', []);

for subj = 1:num_subjects
    
    A_all = [];
    part_coeff_all = [];
    within_mod_z_all = [];
    node_strength_all = [];
    sc_all = [];
    avg_ctrb_disc_all = [];
    mod_ctrb_disc_all = [];
    min_eng_all = [];

    for trial = 1:num_trials
        A_all(:,:,trial) = all_raw_data{subj,trial}.A;
        part_coeff_all(:,trial) = all_raw_data{subj,trial}.part_coeff;
        within_mod_z_all(:,trial) = all_raw_data{subj,trial}.within_mod_z;
        node_strength_all(:,trial) = all_raw_data{subj,trial}.node_strength;
        sc_all(:,trial) = all_raw_data{subj,trial}.sc;
        
        avg_ctrb_disc_all(:,trial) = all_raw_data{subj,trial}.avg_ctrb_disc;
        mod_ctrb_disc_all(:,trial) = all_raw_data{subj,trial}.mod_ctrb_disc;
        min_eng_all(:,trial) = all_raw_data{subj,trial}.min_eng;
    end
    
    data_struct(subj).A = mean(A_all,3);
    data_struct(subj).part_coeff = mean(part_coeff_all,2);
    data_struct(subj).within_mod_z = mean(within_mod_z_all,2);
    data_struct(subj).node_strength = mean(node_strength_all, 2);
    data_struct(subj).sc = mean(sc_all, 2);
    
    data_struct(subj).avg_ctrb_disc = mean(avg_ctrb_disc_all,2);
    data_struct(subj).mod_ctrb_disc = mean(mod_ctrb_disc_all,2); 
    data_struct(subj).min_eng = mean(min_eng_all,2);
end

clearvars -except data_set num_subjects num_trials all_raw_data data_struct

%% Average the Values Across Subjects

part_coeff_all = [];
within_mod_z_all = [];
node_strength_all = [];
sc_all = [];
avg_ctrb_all = [];
mod_ctrb_all = [];
min_eng_all = [];

for subj = 1:num_subjects
    part_coeff_all(:,subj) = data_struct(subj).part_coeff;
    within_mod_z_all(:,subj) = data_struct(subj).within_mod_z;
    node_strength_all(:,subj) = data_struct(subj).node_strength;
    sc_all(:, subj) = data_struct(subj).sc;
    avg_ctrb_all(:,subj) = data_struct(subj).avg_ctrb_disc;
    mod_ctrb_all(:,subj) = data_struct(subj).mod_ctrb_disc;
    min_eng_all(:,subj) = data_struct(subj).min_eng;
end

part_coeff_mean = mean(part_coeff_all,2);
within_mod_z_mean = mean(within_mod_z_all,2);
node_strength_mean = mean(node_strength_all,2);
sc_mean = mean(sc_all,2);

avg_ctrb_mean = mean(avg_ctrb_all,2);
mod_ctrb_mean = mean(mod_ctrb_all,2);
min_eng_mean = mean(min_eng_all,2);

%% Correlations with Degree 

[r_NS_PC, p_NS_PC] = corr(node_strength_mean, part_coeff_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. part_coeff ~ node_strength\n', r_NS_PC, p_NS_PC);
[r_NS_Z, p_NS_Z] = corr(node_strength_mean, within_mod_z_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. within_mod_z ~ node_strength\n', r_NS_Z, p_NS_Z);
[r_NS_SC, p_NS_SC] = corr(node_strength_mean, sc_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. sc ~ node_strength\n', r_NS_SC, p_NS_SC);

%% Minimum Control Energy Models

% Univariate Models

[r_NS_MIN, p_NS_MIN] = corr(node_strength_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ node_strength\n', r_NS_MIN, p_NS_MIN);
[r_PC_MIN, p_PC_MIN] = corr(part_coeff_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ part_coeff\n', r_PC_MIN, p_PC_MIN);
[r_Z_MIN, p_Z_MIN] = corr(within_mod_z_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ within_mod_deg_z\n', r_Z_MIN, p_Z_MIN);
[r_SC_MIN, p_SC_MIN] = corr(sc_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ sc\n', r_SC_MIN, p_SC_MIN);

% Multivariate Models

[r_PC_NS_MIN, p_PC_NS_MIN, residuals] = ... 
    partialcorr_with_resids(part_coeff_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. min_eng ~ part_coeff, controlling for node_strength\n', r_PC_NS_MIN, p_PC_NS_MIN);
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

[r_Z_NS_MIN, p_Z_NS_MIN, residuals] = ... 
    partialcorr_with_resids(within_mod_z_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. min_eng ~ within_mod_deg_z, controlling for node_strength\n', r_Z_NS_MIN, p_Z_NS_MIN);
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

[r_SC_NS_MIN, p_SC_NS_MIN, residuals] = ... 
    partialcorr_with_resids(sc_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. min_eng ~ sc, controlling for node_strength\n', r_SC_NS_MIN, p_SC_NS_MIN);
SC_resid = residuals(:, 1);
MIN_resid = residuals(:, 2);
coeffs = polyfit(SC_resid, MIN_resid, 1);
x = linspace(min(SC_resid), max(SC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SC_resid, MIN_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Subgraph Centrality Residual', 'FontSize', 15);
ylabel('Minimum Control Energy Residual','FontSize',15);
title('SC vs. Minimum Control Energy', ...
    'FontSize', 15);
prettify
hold off

%% Average Controllability Models

% Univariate Models

[r_NS_AVG, p_NS_AVG] = corr(node_strength_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ node_strength\n', r_NS_AVG, p_NS_AVG);
[r_PC_AVG, p_PC_AVG] = corr(part_coeff_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ part_coeff\n', r_PC_AVG, p_PC_AVG);
[r_Z_AVG, p_Z_AVG] = corr(within_mod_z_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ within_mod_deg_z\n', r_Z_AVG, p_Z_AVG);
[r_SC_AVG, p_SC_AVG] = corr(sc_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ sc\n', r_SC_AVG, p_SC_AVG);

% Multivariate Models

[r_PC_NS_AVG, p_PC_NS_AVG, residuals] = ... 
    partialcorr_with_resids(part_coeff_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. avg_ctrb ~ part_coeff, controlling for node_strength\n', r_PC_NS_AVG, p_PC_NS_AVG);
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

[r_Z_NS_AVG, p_Z_NS_AVG, residuals] = ... 
    partialcorr_with_resids(within_mod_z_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. avg_ctrb ~ within_mod_deg_z, controlling for node_strength\n', r_Z_NS_AVG, p_Z_NS_AVG);
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

[r_SC_NS_AVG, p_SC_NS_AVG, residuals] = ... 
    partialcorr_with_resids(sc_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. avg_ctrb ~ sc, controlling for node_strength\n', r_SC_NS_AVG, p_SC_NS_AVG);
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

%% Modal Controllability Models

% Univariate Models

[r_NS_MOD, p_NS_MOD] = corr(node_strength_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ node_strength\n', r_NS_MOD, p_NS_MOD);
[r_PC_MOD, p_PC_MOD] = corr(part_coeff_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ part_coeff\n', r_PC_MOD, p_PC_MOD);
[r_Z_MOD, p_Z_MOD] = corr(within_mod_z_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ within_mod_deg_z\n', r_Z_MOD, p_Z_MOD);
[r_SC_MOD, p_SC_MOD] = corr(sc_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ sc\n', r_SC_MOD, p_SC_MOD);

% Multivariate Models

[r_PC_NS_MOD, p_PC_NS_MOD, residuals] = ... 
    partialcorr_with_resids(part_coeff_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. mod_ctrb ~ part_coeff, controlling for node_strength\n', r_PC_NS_MOD, p_PC_NS_MOD);
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
hold off

[r_Z_NS_MOD, p_Z_NS_MOD, residuals] = ... 
    partialcorr_with_resids(within_mod_z_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. mod_ctrb ~ within_mod_deg_z, controlling for node_strength\n', r_Z_NS_MOD, p_Z_NS_MOD);
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

[r_SC_NS_MOD, p_SC_NS_MOD, residuals] = ... 
    partialcorr_with_resids(sc_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. mod_ctrb ~ sc, controlling for node_strength\n', r_SC_NS_MOD, p_SC_NS_MOD);
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
