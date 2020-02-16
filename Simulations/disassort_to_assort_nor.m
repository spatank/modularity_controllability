%% Set Network Parameters

dist = 'normal';

N = 234; % number of nodes
density = 0.1485; % desired density
K = round(((N^2-N)*density)/2); % number of edges

Ci = zeros(N,1);
Ci(1:234/2) = 1;
Ci((234/2)+1:end) = 2;

%% Set Controllability Parameters

x_0 = 0; % initial state for state transition
x_f = 1; % final state for state transition
T = 1; % time horizon for controllability
nor = 0; % matrix normalization flag for stability
thresh = 1; % threshold for slower and faster modes for modal control

noise_flag = 0;

%% Generate Network Ensembles

iters = 50; % number of networks in each ensemble
vec = linspace(0.001, 0.999, 25); % number of distinct network topologies

Q = zeros(size(vec)); % modularity index
Q_err = zeros(size(vec));
avg_ctrb = zeros(size(vec));
avg_ctrb_err = zeros(size(vec));
mod_ctrb = zeros(size(vec));
mod_ctrb_err = zeros(size(vec));
min_eng = zeros(size(vec));
min_eng_err = zeros(size(vec));

spec_met = zeros(size(vec));

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

all_eig_vals_new = zeros(N, length(vec), iters);
for idx = 1:length(vec)
    for iter = 1:iters
        all_nets(:,:,idx,iter) = (all_nets(:,:,idx,iter)/(max_eig));
        all_eig_vals_new(:,idx,iter) = eig(all_nets(:,:,idx,iter)); % store eigenvalues
    end
end

%% Compute Controllabilities

for idx = 1:length(vec) 
    fprintf('Idx %d\n', idx);
    Q_frac = zeros(1,iters);
    avg_ctrb_frac = zeros(1,iters);
    mod_ctrb_frac = zeros(1,iters);
    min_eng_frac = zeros(1,iters);
    spec_met_frac = zeros(1,iters);
    for iter = 1:iters
        A = all_nets(:,:,idx,iter);
        Q_frac(iter) = modularity_index(A, Ci);
        avg_ctrb_frac(iter) = mean(avg_ctrb_cont(A, T, nor));
        mod_ctrb_frac(iter) = mean(mod_ctrb_cont(A, 1:size(A,1), thresh, nor));
        x0_vec = x_0.*ones(size(A,1),1);
        xf_vec = x_f.*ones(size(A,1),1);
        [~, u, n_err] = min_eng_cont(A, T, eye(size(A)), x0_vec, xf_vec, nor);
        if n_err > 10^-6
            disp('Error Threshold Exceeded')
        end
        del_T = 3/length(u);
        min_eng_frac(iter) = mean(sum(u.^2, 2)*del_T);
        spec_met_frac(iter) = mean(sum(A^2));
    end
    Q(idx) = mean(Q_frac);
    Q_err(idx) = std(Q_frac)/iters;
    avg_ctrb(idx) = mean(avg_ctrb_frac);
    avg_ctrb_err(idx) = std(avg_ctrb_frac)/iters;
    mod_ctrb(idx) = mean(mod_ctrb_frac);
    mod_ctrb_err(idx) = std(mod_ctrb_frac)/iters;
    min_eng(idx) = mean(min_eng_frac);
    min_eng_err(idx) = std(min_eng_frac)/iters;
    spec_met(idx) = mean(spec_met_frac);
end

%% Plots

f = figure('color','w');
errorbar(vec,zscore(avg_ctrb),zscore(avg_ctrb_err),'.k','MarkerSize',25,'LineWidth',0.1);
title('Nor.: Disassort. to Assort.: Avg. Ctrb.', 'FontSize', 15);
xlabel('Fraction of Edges in Core', 'FontSize', 15);
ylabel('Average Controllability', 'FontSize', 15);
prettify

figure;
errorbar(vec,zscore(mod_ctrb),zscore(mod_ctrb_err),'.k','MarkerSize',25,'LineWidth',0.1);
title('Nor.: Disassort. to Assort.: Mod. Ctrb.', 'FontSize', 15);
xlabel('Fraction of Edges in Core', 'FontSize', 15);
ylabel('Modal Controllability', 'FontSize', 15);
prettify

figure;
errorbar(vec,zscore(min_eng), zscore(min_eng_err),'.k','MarkerSize',25,'LineWidth',0.1);
title('Nor.: Disassort. to Assort.: Min. Eng.', 'FontSize', 15);
xlabel('Fraction of Edges in Core', 'FontSize', 15);
ylabel('Control Energy', 'FontSize', 15);
prettify

figure;
errorbar(Q, zscore(avg_ctrb), zscore(avg_ctrb_err),'.k','MarkerSize',25,'LineWidth',0.1);
title('Nor.: Disassort. to Assort.: Avg. Ctrb.', 'FontSize', 15);
xlabel('Modularity Q', 'FontSize', 15);
ylabel('Average Controllability', 'FontSize', 15);
prettify

figure;
errorbar(Q, zscore(mod_ctrb), zscore(mod_ctrb_err),'.k','MarkerSize',25,'LineWidth',0.1);
title('Nor.: Disassort. to Assort.: Mod. Ctrb.', 'FontSize', 15);
xlabel('Modularity Q', 'FontSize', 15);
ylabel('Modal Controllability', 'FontSize', 15);
prettify

figure;
errorbar(Q, zscore(min_eng), zscore(min_eng_err), '.k','MarkerSize',25,'LineWidth',0.1);
title('Nor.: Disassort. to Assort.: Min. Eng.', 'FontSize', 15);
xlabel('Modularity Q', 'FontSize', 15);
ylabel('Control Energy', 'FontSize', 15);
prettify

figure;
plot(vec, Q,'.k','MarkerSize',25,'LineWidth',0.1);
title('Modularity Q and Fraction of Edges in Modules', 'FontSize', 25);
xlabel('Fraction of Edges in Modules', 'FontSize', 15);
ylabel('Modularity Q', 'FontSize', 15);
prettify