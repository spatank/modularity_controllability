function node_level_spectral_metric = node_level_spectral_metric(A)
% node_level_spectral_metric Returns the ratio of the largest to smallest
% eigenvalue-weighted projection of a vector corresponding to a node in the
% network.
% Project the vector corresponding to a node onto the largest eigenmode 
% weighted by the largest eigenvalue. Project the vector corresponding to 
% the node onto the smallest eigenmode. Return the ratio. Perform this 
% operation for all nodes in the network.

% max_eig = eigs(A,1); % pre-normalization
% A = 0.9*(A/(max_eig)) - eye(size(A)); % normalize the adjacency matrix

[V, D] = eig(A);

node_level_spectral_metric = zeros(1, length(A));

for node = 1:size(A, 1)
    node_vec = zeros(length(A), 1);
    node_vec(node) = 1;
    node_level_spectral_metric(node) = ...
        abs(D(1,1)*V(:,1)'*node_vec) - abs(D(end,end)*V(:,end)'*node_vec);
end


end

