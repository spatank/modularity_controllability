function [most_central] = central_partition(Models)
% central_partition Returns the most partitioning of vertices based on
% variation of information (needs WSBM_v1.2)
%   Models: fitted WSBMs
num_partitions = length(Models);
VI_mat = zeros(num_partitions, num_partitions);
for i = 1:num_partitions
    for j = 1:num_partitions
        VI_mat(i, j) = -varInfo(Models{i, 1}.Para.mu, Models{j, 1}.Para.mu);
    end
end
[~, best_model_idx] = max(sum(VI_mat, 2)); % index of most central model
[~, most_central] = max(Models{best_model_idx, 1}.Para.mu);
end

