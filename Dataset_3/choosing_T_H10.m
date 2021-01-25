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
T_rng = 1:1:50; % time horizon for controllability
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

r_NS_PC = zeros(size(T_rng));
p_NS_PC = zeros(size(T_rng));
r_NS_Z = zeros(size(T_rng));
p_NS_Z = zeros(size(T_rng));
r_NS_SPEC = zeros(size(T_rng));
p_NS_SPEC = zeros(size(T_rng));

r_NS_MIN = zeros(size(T_rng));
p_NS_MIN = zeros(size(T_rng));
r_PC_MIN = zeros(size(T_rng));
p_PC_MIN = zeros(size(T_rng));
r_Z_MIN = zeros(size(T_rng));
p_Z_MIN = zeros(size(T_rng));
r_SC_MIN = zeros(size(T_rng));
p_SC_MIN = zeros(size(T_rng));
r_PC_NS_MIN = zeros(size(T_rng)); 
p_PC_NS_MIN = zeros(size(T_rng));
r_Z_NS_MIN = zeros(size(T_rng)); 
p_Z_NS_MIN = zeros(size(T_rng));
r_SC_NS_MIN = zeros(size(T_rng));
p_SC_NS_MIN = zeros(size(T_rng));

r_NS_AVG = zeros(size(T_rng));
p_NS_AVG = zeros(size(T_rng));
r_PC_AVG = zeros(size(T_rng));
p_PC_AVG = zeros(size(T_rng));
r_Z_AVG = zeros(size(T_rng));
p_Z_AVG = zeros(size(T_rng));
r_SC_AVG = zeros(size(T_rng));
p_SC_AVG = zeros(size(T_rng));
r_PC_NS_AVG = zeros(size(T_rng)); 
p_PC_NS_AVG = zeros(size(T_rng));
r_Z_NS_AVG = zeros(size(T_rng)); 
p_Z_NS_AVG = zeros(size(T_rng));
r_SC_NS_AVG = zeros(size(T_rng));
p_SC_NS_AVG = zeros(size(T_rng));

r_NS_MOD = zeros(size(T_rng));
p_NS_MOD = zeros(size(T_rng));
r_PC_MOD = zeros(size(T_rng));
p_PC_MOD = zeros(size(T_rng));
r_Z_MOD = zeros(size(T_rng));
p_Z_MOD = zeros(size(T_rng));
r_SC_MOD = zeros(size(T_rng));
p_SC_MOD = zeros(size(T_rng));
r_PC_NS_MOD = zeros(size(T_rng)); 
p_PC_NS_MOD = zeros(size(T_rng));
r_Z_NS_MOD = zeros(size(T_rng)); 
p_Z_NS_MOD = zeros(size(T_rng));
r_SC_NS_MOD = zeros(size(T_rng));
p_SC_NS_MOD = zeros(size(T_rng));

for T_idx = 1:length(T_rng)
    
    all_raw_data = cell(1, num_subjects);
    
    disp(['T = ', num2str(T_rng(T_idx))]);
    
    T_avg = T_rng(T_idx); % time horizon for average controllability and minimum control energy
    T_eng = T_rng(T_idx); % time horizon for minimum control energy

    for subj = 1:num_subjects

        subj_struct = struct('A', [], 'M', [], ...
        'part_coeff', [], 'within_mod_z', [], ...
        'node_strength', [], 'sc', [], ...
        'avg_ctrb_disc',[], 'mod_ctrb_disc',[], ...
        'min_eng', []);

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

        subj_struct.sc = diag(expm(A)); % subgraph centrality

        subj_struct.avg_ctrb_disc = avg_ctrb_disc(A, T_avg, nor);
        subj_struct.mod_ctrb_disc = mod_ctrb_disc(A, 1:size(A,1), thresh, nor);

        subj_struct.min_eng = min_eng_0_1_node(A, T_eng, 'disc');

        all_raw_data{subj} = subj_struct;
    end

    % Average the Values Across the Subjects

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
        sc_all(:, subj) = all_raw_data{1, subj}.sc;
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

    % Correlations with Degree 

    [r_NS_PC(T_idx), p_NS_PC(T_idx)] = corr(node_strength_mean, part_coeff_mean, 'type', 'Spearman');
    [r_NS_Z(T_idx), p_NS_Z(T_idx)] = corr(node_strength_mean, within_mod_z_mean, 'type', 'Spearman');
    [r_NS_SPEC(T_idx), p_NS_SPEC(T_idx)] = corr(node_strength_mean, sc_mean, 'type', 'Spearman');

    % Minimum Control Energy Models

    % Univariate Models
    [r_NS_MIN(T_idx), p_NS_MIN(T_idx)] = corr(node_strength_mean, min_eng_mean, 'type', 'Spearman');
    [r_PC_MIN(T_idx), p_PC_MIN(T_idx)] = corr(part_coeff_mean, min_eng_mean, 'type', 'Spearman');
    [r_Z_MIN(T_idx), p_Z_MIN(T_idx)] = corr(within_mod_z_mean, min_eng_mean, 'type', 'Spearman');
    [r_SC_MIN(T_idx), p_SC_MIN(T_idx)] = corr(sc_mean, min_eng_mean, 'type', 'Spearman');

    % Multivariate Models
    [r_PC_NS_MIN(T_idx), p_PC_NS_MIN(T_idx), residuals] = ... 
        partialcorr_with_resids(part_coeff_mean, min_eng_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');
    [r_Z_NS_MIN(T_idx), p_Z_NS_MIN(T_idx), residuals] = ... 
        partialcorr_with_resids(within_mod_z_mean, min_eng_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');
    [r_SC_NS_MIN(T_idx), p_SC_NS_MIN(T_idx), residuals] = ... 
        partialcorr_with_resids(sc_mean, min_eng_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');

    % Average Controllability Models

    % Univariate Models

    [r_NS_AVG(T_idx), p_NS_AVG(T_idx)] = corr(node_strength_mean, avg_ctrb_mean, 'type', 'Spearman');
    [r_PC_AVG(T_idx), p_PC_AVG(T_idx)] = corr(part_coeff_mean, avg_ctrb_mean, 'type', 'Spearman');
    [r_Z_AVG(T_idx), p_Z_AVG(T_idx)] = corr(within_mod_z_mean, avg_ctrb_mean, 'type', 'Spearman');
    [r_SC_AVG(T_idx), p_SC_AVG(T_idx)] = corr(sc_mean, avg_ctrb_mean, 'type', 'Spearman');

    % Multivariate Models

    [r_PC_NS_AVG(T_idx), p_PC_NS_AVG(T_idx), residuals] = ... 
        partialcorr_with_resids(part_coeff_mean, avg_ctrb_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');
    [r_Z_NS_AVG(T_idx), p_Z_NS_AVG(T_idx), residuals] = ... 
        partialcorr_with_resids(within_mod_z_mean, avg_ctrb_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');
    [r_SC_NS_AVG(T_idx), p_SC_NS_AVG(T_idx), residuals] = ... 
        partialcorr_with_resids(sc_mean, avg_ctrb_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');

    % Modal Controllability Models

    % Univariate Models

    [r_NS_MOD(T_idx), p_NS_MOD(T_idx)] = corr(node_strength_mean, mod_ctrb_mean, 'type', 'Spearman');
    [r_PC_MOD(T_idx), p_PC_MOD(T_idx)] = corr(part_coeff_mean, mod_ctrb_mean, 'type', 'Spearman');
    [r_Z_MOD(T_idx), p_Z_MOD(T_idx)] = corr(within_mod_z_mean, mod_ctrb_mean, 'type', 'Spearman');
    [r_SC_MOD(T_idx), p_SC_MOD(T_idx)] = corr(sc_mean, mod_ctrb_mean, 'type', 'Spearman');

    % Multivariate Models

    [r_PC_NS_MOD(T_idx), p_PC_NS_MOD(T_idx), residuals] = ... 
        partialcorr_with_resids(part_coeff_mean, mod_ctrb_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');
    [r_Z_NS_MOD(T_idx), p_Z_NS_MOD(T_idx), residuals] = ... 
        partialcorr_with_resids(within_mod_z_mean, mod_ctrb_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');
    [r_SC_NS_MOD(T_idx), p_SC_NS_MOD(T_idx), residuals] = ... 
        partialcorr_with_resids(sc_mean, mod_ctrb_mean, node_strength_mean, ...
        'type', 'Spearman', 'rows', 'complete');

end

%% Make Figures

path_1 = '/Users/sppatankar/Desktop/Projects/Modularity/Re-submission/Choosing_T/';

f_1 = figure('color', 'w');
hold on
plot(T_rng, r_SC_NS_MIN, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 25, 'LineWidth', 0.1);
xlabel('T', 'FontSize', 15);
ylabel('r(SC, MIN) corrected for NS', 'FontSize', 15);
title('r(SC, MIN) vs. T', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('r(SC, MIN)_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

f_2 = figure('color', 'w');
hold on
plot(T_rng, r_SC_NS_AVG, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 25, 'LineWidth', 0.1);
xlabel('T', 'FontSize', 15);
ylabel('r(SC, AVG) corrected for NS', 'FontSize', 15);
title('r(SC, AVG) vs. T', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('r(SC, AVG)_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

% f_3 = figure('color', 'w');
% hold on
% plot(T_rng, r_SC_NS_MOD, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 25, 'LineWidth', 0.1);
% xlabel('T', 'FontSize', 15);
% ylabel('r(SC, MOD) corrected for NS', 'FontSize', 15);
% title('r(SC, MOD) vs. T', ...
%     'FontSize', 15);
% prettify
% hold off
% path_2 = strcat('r(SC, MOD)_', data_set);
% saveas(gcf, fullfile(path_1, path_2), 'epsc')

f_4 = figure('color', 'w');
hold on
plot(T_rng, p_SC_NS_MIN, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 25, 'LineWidth', 0.1);
xlabel('T', 'FontSize', 15);
ylabel('p(SC, MIN) corrected for NS', 'FontSize', 15);
title('p(SC, MIN) vs. T', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('p(SC, MIN)_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

f_5 = figure('color', 'w');
hold on
plot(T_rng, p_SC_NS_AVG, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 25, 'LineWidth', 0.1);
xlabel('T', 'FontSize', 15);
ylabel('p(SC, AVG) corrected for NS', 'FontSize', 15);
title('p(SC, AVG) vs. T', ...
    'FontSize', 15);
prettify
hold off
path_2 = strcat('p(SC, AVG)_', data_set);
saveas(gcf, fullfile(path_1, path_2), 'epsc')

% f_6 = figure('color', 'w');
% hold on
% plot(T_rng, p_SC_NS_MOD, '.', 'MarkerFaceColor', 'k', 'MarkerSize', 25, 'LineWidth', 0.1);
% xlabel('T', 'FontSize', 15);
% ylabel('p(SC, MOD) corrected for NS', 'FontSize', 15);
% title('p(SC, MOD) vs. T', ...
%     'FontSize', 15);
% prettify
% hold off
% path_2 = strcat('p(SC, MOD)_', data_set);
% saveas(gcf, fullfile(path_1, path_2), 'epsc')


