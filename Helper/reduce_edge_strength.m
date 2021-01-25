function A_new = reduce_edge_strength(A, percent_edges, reduce_by)
% reduce_edge_strength Reduces the edge strength of a pre-specified
% percentage of total edges by a pre-specified percent amount.
%   A: weighted adjacency matrix
%   percent_edges: percentage of edges to reduce
%   reduce_by: percent to reduce by

A(A == 0) = NaN;
edge_list = Adj2Edg(A);
x_percent = round((percent_edges/100) * length(edge_list)); % how many edges to reduce
selected_edges = randsample(length(edge_list), x_percent);

for idx = 1:length(selected_edges)
    selected_edge = selected_edges(idx);
    edge_list(selected_edge,3) = edge_list(selected_edge,3) - ...
        ((reduce_by/100) * edge_list(selected_edge,3));
end

A_new = Edg2Adj(edge_list);

A_new(isnan(A_new)) = 0;

end

