function [ phi_p, phi_t ] = mod_ctrb_disc(A, i, thresh, nor)
% adapted from Jeni's code
% Calculate persistent and transient modal controllability of node i. 
% Intuitively, nodes with high persistent controllability will result in
% large perturbations to slow modes of the system when stimulated, and
% nodes with high transient controllability will result in large
% perturbations to fast modes of the system when stimulated
% Inputs:  
% A         The NxN structural (NOT FUNCTIONAL) network adjacency matrix
% i         This index of the node you wish to calculate the
%           controllability of
% thresh    The thresh hold for calculating persistent and transient modal
%           controllability. For example, if thresh is .1, this will
%           use the 10% slowest modes for persistent controllability, and
%           the 10% fastest modes for transient controllability.
% nor       normalization constant (default: 1+eig_largest) (SP)

% Outputs:
% phi_p     Persistent controllability of the node (should be between 0 and
%           1)
% phi_t     Transient controllability (should be between 0 and 1)

% @author JStiso

% Normalize
% if nor==1
%     A = A./(1+eigs(A,1));
% end

if nargin ~= 4
    nor = (1+eigs(A,1));
    A = A./nor;
else
    A = A./nor;
end

[V, D] = eig(A);
lambda = eig(A);

% get only the X% largest values
lambda = abs(lambda);
[modes, idx] = sort(lambda);
cutoff = round(numel(modes)*thresh);
transient_modes = modes(1:cutoff);
persistent_modes = modes((end-cutoff+1):end);
% get controllability
phi_p = sum((1 - (persistent_modes).^2).*(V(i,idx((numel(idx)-cutoff+1):end)).^2)');
phi_t = sum((1 - (transient_modes).^2).*(V(i,idx(1:cutoff)).^2)');

end
