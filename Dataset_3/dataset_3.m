clear; close all;

%% Generate Data Structure

% input parameters for controllability metrics
T = 1; % time horizon for average controllability and minimum control energy
x_0 = 0; % initial state for state transition
x_f = 1; % final state for state transition
nor = 0; % matrix normalization flag for stability
thresh = 1; % threshold for slower and faster modes for modal control

path_1 = '/Volumes/Elements/Modularity/Human_10_x2/WSBM_Results/'; % LN
% path_1 = '/Volumes/Elements/Modularity/Human_10/WSBM_Results/'; % Normal

num_subjects = 10;

all_raw_data = cell(1,num_subjects);

for subj = 1:num_subjects

    subj_struct = struct('A', [], 'M', [], ...
    'part_coeff', [], 'within_mod_z', [], ...
    'node_strength', [], 'communicability', [], ...
    'spec_met', [], ...
    'avg_ctrb_cont',[], 'mod_ctrb_cont',[], ...
    'min_eng', []);

    fprintf('Subj = %d.\n', subj);
    path_2 = sprintf('Data_Copy/Subj_%d.mat', subj);
    load(fullfile(path_1, path_2)); % loads 'connectivity' and 'name'
    
    A = connectivity;
    
    % A = 0.9*(A/(eigs(A,1)));
    A = (A/(eigs(A,1)));
    
    subj_struct.A = A;
    
    k = 10; % number of communities; with LogNormal prior
    % k = 11; % number of communities; with Normal prior
    
    path_3 = sprintf('k_%d/subj_%d_k_%d_human_10_LN.mat', k, subj, k);
    % path_3 = sprintf('k_%d/subj_%d_k_%d_human_10.mat', k, subj, k);
    load(fullfile(path_1, path_3)); % loads 'Best_Model', 'Models', and 'Scores'

    Ci = zeros(size(A, 1), length(Models)); % multiple fitted WSBMs 
    for idx = 1:length(Models)
        model = Models(idx);
        labels = assign_communities(model{1, 1});
        Ci(:, idx) = labels;
    end
    
    M = central_partition(Ci);
    % M = M(1:233);
    subj_struct.M = M;

    subj_struct.part_coeff = participation_coef(A, M);
    subj_struct.within_mod_z = module_degree_zscore(A, M);
    subj_struct.node_strength = strengths_und(A);
    
    subj_struct.communicability = sum(expm(A)); % communicability
    subj_struct.spec_met = node_level_spectral_metric(A);
    
    subj_struct.avg_ctrb_cont = avg_ctrb_cont(A, T, nor);
    subj_struct.mod_ctrb_cont = mod_ctrb_cont(A, 1:size(A, 1), thresh, nor);
    
    % compute minimum control energy
    x0_vec = x_0.*ones(size(A, 1), 1);
    xf_vec = x_f.*ones(size(A, 1), 1);
    [~, u, n_err] = min_eng_cont(A, T, eye(size(A)), x0_vec, xf_vec, nor);
    if n_err > 10^-6
        disp('Error: Threshold Exceeded')
    end
    del_T = 3/length(u);
    subj_struct.min_eng = sum(u.^2, 2)*del_T;
   
    
    all_raw_data{subj} = subj_struct;
end

%% Average the Values Across Subjects

part_coeff_all = [];
within_mod_z_all = [];
node_strength_all = [];
communicability_all = [];
spec_met_all = [];
avg_ctrb_all = [];
mod_ctrb_all = [];
min_eng_all = [];

for subj = 1:num_subjects
    part_coeff_all(:, subj) = all_raw_data{1, subj}.part_coeff;
    within_mod_z_all(:, subj) = all_raw_data{1, subj}.within_mod_z;
    node_strength_all(:, subj) = all_raw_data{1, subj}.node_strength;
    communicability_all(:, subj) = all_raw_data{1, subj}.communicability;
    spec_met_all(:, subj) = all_raw_data{1, subj}.spec_met;
    avg_ctrb_all(:, subj) = all_raw_data{1, subj}.avg_ctrb_cont;
    mod_ctrb_all(:, subj) = all_raw_data{1, subj}.mod_ctrb_cont;
    min_eng_all(:, subj) = all_raw_data{1, subj}.min_eng;
end

part_coeff_mean = mean(part_coeff_all, 2);
within_mod_z_mean = mean(within_mod_z_all, 2);
node_strength_mean = mean(node_strength_all, 2);
communicability_mean = mean(communicability_all, 2);
spec_met_mean = mean(spec_met_all, 2);

avg_ctrb_mean = mean(avg_ctrb_all, 2);
mod_ctrb_mean = mean(mod_ctrb_all, 2);
min_eng_mean = mean(min_eng_all, 2);


%% Correlations with Degree 

[r_NS_PC, p_NS_PC] = corr(node_strength_mean, part_coeff_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. part_coeff ~ node_strength\n', r_NS_PC, p_NS_PC);
[r_NS_Z, p_NS_Z] = corr(node_strength_mean, within_mod_z_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. within_mod_z ~ node_strength\n', r_NS_Z, p_NS_Z);
% [r_NS_COM, p_NS_COM] = corr(node_strength_mean, communicability_mean, 'type', 'Spearman');
% fprintf('r = %f, p = %f. communicability ~ node_strength\n', r_NS_COM, p_NS_COM);
[r_NS_SPEC, p_NS_SPEC] = corr(node_strength_mean, spec_met_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. communicability ~ node_strength\n', r_NS_SPEC, p_NS_SPEC);

fprintf('\n');

%% Minimum Control Energy Models

% Univariate Models

[r_NS_MIN, p_NS_MIN] = corr(node_strength_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ node_strength\n', r_NS_MIN, p_NS_MIN);
[r_PC_MIN, p_PC_MIN] = corr(part_coeff_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ part_coeff\n', r_PC_MIN, p_PC_MIN);
[r_Z_MIN, p_Z_MIN] = corr(within_mod_z_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ within_mod_deg_z\n', r_Z_MIN, p_Z_MIN);
% [r_COM_MIN, p_COM_MIN] = corr(communicability_mean, min_eng_mean, 'type', 'Spearman');
% fprintf('r = %f, p = %f. min_eng ~ communicability\n', r_COM_MIN, p_COM_MIN);
[r_SPEC_MIN, p_SPEC_MIN] = corr(spec_met_mean, min_eng_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. min_eng ~ spec_met\n', r_SPEC_MIN, p_SPEC_MIN);

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

% [r_COM_NS_MIN, p_COM_NS_MIN, residuals] = ... 
%     partialcorr_with_resids(communicability_mean, min_eng_mean, node_strength_mean, ...
%     'type', 'Spearman', 'rows', 'complete');
% fprintf('r = %f, p = %f. min_eng ~ communicability, controlling for node_strength\n', r_COM_NS_MIN, p_COM_NS_MIN);
% COM_resid = residuals(:, 1);
% MIN_resid = residuals(:, 2);
% coeffs = polyfit(COM_resid, MIN_resid, 1);
% x = linspace(min(COM_resid), max(COM_resid), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on
% plot(COM_resid, MIN_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
% plot(x, y, 'LineWidth', 2);
% xlabel('Communicability Residual','FontSize',15);
% ylabel('Minimum Control Energy Residual','FontSize',15);
% title('Communicability vs. Minimum Control Energy', ...
%     'FontSize', 15);
% prettify
% hold off

[r_SPEC_NS_MIN, p_SPEC_NS_MIN, residuals] = ... 
    partialcorr_with_resids(spec_met_mean, min_eng_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. min_eng ~ spec_met, controlling for node_strength\n', r_SPEC_NS_MIN, p_SPEC_NS_MIN);
SPEC_resid = residuals(:, 1);
MIN_resid = residuals(:, 2);
coeffs = polyfit(SPEC_resid, MIN_resid, 1);
x = linspace(min(SPEC_resid), max(SPEC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SPEC_resid, MIN_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Spec. Met. Residual','FontSize',15);
ylabel('Minimum Control Energy Residual','FontSize',15);
title('Spec. Met. vs. Minimum Control Energy', ...
    'FontSize', 15);
prettify
hold off

fprintf('\n');

%% Average Controllability Models

% Univariate Models

[r_NS_AVG, p_NS_AVG] = corr(node_strength_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ node_strength\n', r_NS_AVG, p_NS_AVG);
[r_PC_AVG, p_PC_AVG] = corr(part_coeff_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ part_coeff\n', r_PC_AVG, p_PC_AVG);
[r_Z_AVG, p_Z_AVG] = corr(within_mod_z_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ within_mod_deg_z\n', r_Z_AVG, p_Z_AVG);
% [r_COM_AVG, p_COM_AVG] = corr(communicability_mean, avg_ctrb_mean, 'type', 'Spearman');
% fprintf('r = %f, p = %f. avg_ctrb ~ communicability\n', r_COM_AVG, p_COM_AVG);
[r_SPEC_AVG, p_SPEC_AVG] = corr(spec_met_mean, avg_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. avg_ctrb ~ spec_met\n', r_SPEC_AVG, p_SPEC_AVG);

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

% [r_COM_NS_AVG, p_COM_NS_AVG, residuals] = ... 
%     partialcorr_with_resids(communicability_mean, avg_ctrb_mean, node_strength_mean, ...
%     'type', 'Spearman', 'rows', 'complete');
% fprintf('r = %f, p = %f. avg_ctrb_mean ~ communicability, controlling for node_strength\n', r_COM_NS_AVG, p_COM_NS_AVG);
% COM_resid = residuals(:, 1);
% AVG_resid = residuals(:, 2);
% coeffs = polyfit(COM_resid, AVG_resid, 1);
% x = linspace(min(COM_resid), max(COM_resid), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on
% plot(COM_resid, AVG_resid,'.','MarkerFaceColor','k','MarkerSize',15,'LineWidth',2);
% plot(x, y, 'LineWidth', 2);
% xlabel('Communicability Residual', 'FontSize', 15);
% ylabel('Average Controllability Residual', 'FontSize', 15);
% title('Communicability vs. Average Controllability', ...
%     'FontSize', 15);
% prettify
% hold off

[r_SPEC_NS_AVG, p_SPEC_NS_AVG, residuals] = ... 
    partialcorr_with_resids(spec_met_mean, avg_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. avg_ctrb ~ spec_met, controlling for node_strength\n', r_SPEC_NS_AVG, p_SPEC_NS_AVG);
SPEC_resid = residuals(:, 1);
AVG_resid = residuals(:, 2);
coeffs = polyfit(SPEC_resid, AVG_resid, 1);
x = linspace(min(SPEC_resid), max(SPEC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SPEC_resid, AVG_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Spec. Met. Residual', 'FontSize', 15);
ylabel('Average Controllability Residual', 'FontSize', 15);
title('Spec. Met. vs. Average Controllability', ...
    'FontSize', 15);
prettify
hold off

fprintf('\n');

%% Modal Controllability Models

% Univariate Models

[r_NS_MOD, p_NS_MOD] = corr(node_strength_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ node_strength\n', r_NS_MOD, p_NS_MOD);
[r_PC_MOD, p_PC_MOD] = corr(part_coeff_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ part_coeff\n', r_PC_MOD, p_PC_MOD);
[r_Z_MOD, p_Z_MOD] = corr(within_mod_z_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ within_mod_deg_z\n', r_Z_MOD, p_Z_MOD);
% [r_COM_MOD, p_COM_MOD] = corr(communicability_mean, mod_ctrb_mean, 'type', 'Spearman');
% fprintf('r = %f, p = %f. mod_ctrb ~ communicability\n', r_COM_MOD, p_COM_MOD);
[r_SPEC_MOD, p_SPEC_MOD] = corr(spec_met_mean, mod_ctrb_mean, 'type', 'Spearman');
fprintf('r = %f, p = %f. mod_ctrb ~ spec_met\n', r_SPEC_MOD, p_SPEC_MOD);

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

% [r_COM_NS_MOD, p_COM_NS_MOD, residuals] = ... 
%     partialcorr_with_resids(communicability_mean, mod_ctrb_mean, node_strength_mean, ...
%     'type', 'Spearman', 'rows', 'complete');
% fprintf('r = %f, p = %f. mod_ctrb ~ communicability, controlling for node_strength\n', r_COM_NS_MOD, p_COM_NS_MOD);
% COM_resid = residuals(:, 1);
% MOD_resid = residuals(:, 2);
% coeffs = polyfit(COM_resid, MOD_resid, 1);
% x = linspace(min(COM_resid), max(COM_resid), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on
% plot(COM_resid, MOD_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
% plot(x, y, 'LineWidth', 2);
% xlabel('Communicability Residual', 'FontSize', 15);
% ylabel('Modal Controllability Residual', 'FontSize', 15);
% title('Communicability vs. Modal Controllability', ...
%     'FontSize', 15);
% prettify
% hold off

[r_SPEC_NS_MOD, p_SPEC_NS_MOD, residuals] = ... 
    partialcorr_with_resids(spec_met_mean, mod_ctrb_mean, node_strength_mean, ...
    'type', 'Spearman', 'rows', 'complete');
fprintf('r = %f, p = %f. mod_ctrb ~ spec_met, controlling for node_strength\n', r_SPEC_NS_MOD, p_SPEC_NS_MOD);
SPEC_resid = residuals(:, 1);
MOD_resid = residuals(:, 2);
coeffs = polyfit(SPEC_resid, MOD_resid, 1);
x = linspace(min(SPEC_resid), max(SPEC_resid), 1000);
y = polyval(coeffs, x);
f = figure('color', 'w');
hold on
plot(SPEC_resid, MOD_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
plot(x, y, 'LineWidth', 2);
xlabel('Spec. Met. Residual', 'FontSize', 15);
ylabel('Modal Controllability Residual', 'FontSize', 15);
title('Spec. Met. vs. Modal Controllability', ...
    'FontSize', 15);
prettify
hold off

% %% Communicability 
% 
% [r_COM_NS_PC, p_COM_NS_PC, residuals] = ... 
%     partialcorr_with_resids(communicability_mean, part_coeff_mean, node_strength_mean, ...
%     'type', 'Spearman', 'rows', 'complete');
% fprintf('r = %f, p = %f. part_coeff ~ communicability, controlling for node_strength\n', r_COM_NS_PC, p_COM_NS_PC);
% COM_resid = residuals(:, 1);
% PC_resid = residuals(:, 2);
% coeffs = polyfit(COM_resid, PC_resid, 1);
% x = linspace(min(COM_resid), max(COM_resid), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on
% plot(COM_resid, PC_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
% plot(x, y, 'LineWidth', 2);
% xlabel('Communicability Residual', 'FontSize', 15);
% ylabel('Participation Coefficient Residual','FontSize',15);
% title('Communicability vs. Participation Coefficient', ...
%     'FontSize', 15);
% prettify
% hold off
% 
% [r_COM_NS_Z, p_COM_NS_Z, residuals] = ... 
%     partialcorr_with_resids(communicability_mean, within_mod_z_mean, node_strength_mean, ...
%     'type', 'Spearman', 'rows', 'complete');
% fprintf('r = %f, p = %f. within_mod_z ~ communicability, controlling for node_strength\n', r_COM_NS_Z, p_COM_NS_Z);
% COM_resid = residuals(:, 1);
% Z_resid = residuals(:, 2);
% coeffs = polyfit(COM_resid, Z_resid, 1);
% x = linspace(min(COM_resid), max(COM_resid), 1000);
% y = polyval(coeffs, x);
% f = figure('color', 'w');
% hold on
% plot(COM_resid, Z_resid, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 15, 'LineWidth', 2);
% plot(x, y, 'LineWidth', 2);
% xlabel('Communicability Residual','FontSize',15);
% ylabel('Within-Module Degree Z-Score Residual', 'FontSize', 15);
% title('Communicability vs. Within-Module Degree Z-Score', ...
%     'FontSize', 15);
% prettify
% hold off
