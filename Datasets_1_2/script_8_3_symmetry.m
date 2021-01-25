clc; clear; close all;

data_set = 'd_2_1';

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
T_avg = 4; % time horizon for average controllability 
T_eng = 4; % time horizon for minimum control energy
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

all_raw_data = cell(num_subjects,num_trials);

for subj = 1:num_subjects

    subj_struct = struct('A',[],...
    'node_strength',[],'sc',[],...
    'avg_ctrb_disc',[],'min_eng',[]);
    
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
        
        A_raw = reduce_edge_strength(A_raw, 50, 100); % reduce
        
        A = (A_raw/(eigs(A_raw,1)));
        
        subj_struct.A = A;

        subj_struct.node_strength = strengths_und(A);
 
        subj_struct.sc = diag(expm(A)); % communicability
        
        subj_struct.avg_ctrb_disc = avg_ctrb_disc(A, T_avg, nor);
        subj_struct.min_eng = min_eng_0_1_node(A, T_eng, 'disc');

        all_raw_data{subj,trial} = subj_struct;
    end
end


%% Average the Values Across the Three Scanning Sessions

data_struct = struct('A', [], ...
    'node_strength', [],'sc', [], ...
    'avg_ctrb_disc', [], 'min_eng', []);

for subj = 1:num_subjects
    
    A_all = [];
    node_strength_all = [];
    sc_all = [];
    avg_ctrb_disc_all = [];
    min_eng_all = [];

    for trial = 1:num_trials
        A_all(:,:,trial) = all_raw_data{subj,trial}.A;
        node_strength_all(:,trial) = all_raw_data{subj,trial}.node_strength;
        sc_all(:,trial) = all_raw_data{subj,trial}.sc;
        
        avg_ctrb_disc_all(:,trial) = all_raw_data{subj,trial}.avg_ctrb_disc;
        min_eng_all(:,trial) = all_raw_data{subj,trial}.min_eng;
    end
    
    data_struct(subj).A = mean(A_all,3);
    data_struct(subj).node_strength = mean(node_strength_all, 2);
    data_struct(subj).sc = mean(sc_all, 2);
    
    data_struct(subj).avg_ctrb_disc = mean(avg_ctrb_disc_all,2);
    data_struct(subj).min_eng = mean(min_eng_all,2);
end

clearvars -except data_set num_subjects num_trials all_raw_data data_struct

%% Average the Values Across Subjects

node_strength_all = [];
sc_all = [];
avg_ctrb_all = [];
min_eng_all = [];

for subj = 1:num_subjects
    node_strength_all(:,subj) = data_struct(subj).node_strength;
    sc_all(:, subj) = data_struct(subj).sc;
    avg_ctrb_all(:,subj) = data_struct(subj).avg_ctrb_disc;
    min_eng_all(:,subj) = data_struct(subj).min_eng;
end

node_strength_mean = mean(node_strength_all,2);
sc_mean = mean(sc_all,2);

avg_ctrb_mean = mean(avg_ctrb_all,2);
min_eng_mean = mean(min_eng_all,2);

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
xlabel('Subgraph Centrality Residual', 'FontSize', 15);
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