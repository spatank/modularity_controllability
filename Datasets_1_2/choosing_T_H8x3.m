clc; clear; close all;

data_set = 'd_2_2';

if strcmp(data_set, 'd_1')
    load('ncomm_8_3_qa_gfa.mat'); % for QA weighted streamlines
end
if strcmp(data_set, 'd_2_1') || strcmp(data_set, 'd_2_2')
    load('data8x3fullscale.mat'); % for streamlines
end

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

if strcmp(data_set, 'd_1') % QA with normal prior
    all_As = qaNetworks(:,:,3); % Scale 125: 234 regions
    path_1 = ['/Volumes/My Passport/Modularity_2/Human_8_3_x2/' ...
        'Data_Scripts/WSBM_Results']; 
    k = 12; % number of communities for QA  with Gaussian prior
end
if strcmp(data_set, 'd_2_1') % streamlines with log-normal prior
    all_As = A{1,3}; % Scale 125: 234 regions
    path_1 = ['/Volumes/My Passport/Modularity_2/Human_8_3_x2/' ...
        'Data_Scripts_Streamlines_LogNormal/WSBM_Results_Streamlines_LN'];
    k = 14; % number of communities for streamlines with LN prior
end
if strcmp(data_set, 'd_2_2') % streamlines with Gaussian prior
    all_As = A{1,3}; % Scale 125: 234 regions
    path_1 = ['/Volumes/My Passport/Modularity_2/Human_8_3_x2/' ...
        'Data_Scripts_Streamlines_Normal/WSBM_Results_Streamlines_Normal']; 
    k = 12; % number of communities for streamlines with Gaussian prior
end

num_subjects = size(all_As,1);
num_trials = size(all_As,2);

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
    
    all_raw_data = cell(num_subjects,num_trials);
    
    disp(['T = ', num2str(T_rng(T_idx))]);
    
    T_avg = T_rng(T_idx); % time horizon for average controllability and minimum control energy
    T_eng = T_rng(T_idx); % time horizon for minimum control energy

    for subj = 1:num_subjects

        subj_struct = struct('A',[],'M',[],...
        'part_coeff',[],'within_mod_z',[],...
        'node_strength',[],'sc',[],...
        'avg_ctrb_disc',[],'mod_ctrb_disc',[],...
        'min_eng',[]);

        for trial = 1:num_trials

            path_2 = sprintf('Trial %d', trial);
            
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

            A = (A_raw/(eigs(A_raw,1)));

            subj_struct.A = A;

            load(fullfile(path_1, path_2, path_3)); % loads 'Best_Model', 'Models', and 'Scores'

            Ci = zeros(size(A, 1), length(Models)); % multiple fitted WSBMs 
            for idx = 1:length(Models)
                model = Models(idx);
                labels = assign_communities(model{1,1});
                Ci(:, idx) = labels;
            end

            M = central_partition(Ci);
            subj_struct.M = M;

            subj_struct.part_coeff = participation_coef(A, M);
            subj_struct.within_mod_z = module_degree_zscore(A, M);
            subj_struct.node_strength = strengths_und(A);

            % subj_struct.sc = node_level_spectral_metric(A); 
            % subj_struct.sc = diag(A^2); 
            subj_struct.sc = diag(expm(A)); % subgraph centrality

            subj_struct.avg_ctrb_disc = avg_ctrb_disc(A, T_avg, nor);
            subj_struct.mod_ctrb_disc = mod_ctrb_disc(A, 1:size(A,1), thresh, nor);

            subj_struct.min_eng = min_eng_0_1_node(A, T_eng, 'disc');


            all_raw_data{subj,trial} = subj_struct;
        end
    end

    % Average the Values Across the Three Scanning Sessions

    data_struct = struct('A', [], ...
        'part_coeff', [], 'within_mod_z', [], ...
        'node_strength', [], 'sc', [], ...
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

    % Average the Values Across Subjects

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


