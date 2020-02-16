function [ phi_p, phi_t ] = mod_ctrb_cont(A, i, thresh, nor)
% Calculate persistent and transient modal controllability of node i. 
% Intuitively, nodes with high persistent controllability will result in
% large perturbations to slow modes of the system when stimulated, and
% nodes with high transient controllability will result in large
% perturbations to fast modes of the system when stimulated
% Inputs:  
% A         The NxN structural (NOT FUNCTIONAL) network adjacency matrix
% V         The NxN matrix conatining eigenvectors of your system, where N
%           is the number of nodes
% lambda    The Nx1 vector of eigenvalues for your system
% i         This index of the node you wish to calculate the
%           controllability of
% dt        The time step of the system
% thresh    The thresh hold for calculating persistent and transient modal
%           controllability. For example, if thresh is .1, this will
%           use the 10% slowest modes for persistent controllability, and
%           the 10% fastest modes for transient controllability.
% nor       normalization constant (SP)

% Outputs:
% phi_p     Persistent controllability of the node (should be between 0 and
%           1)
% phi_t     Transient controllability (should be between 0 and 1)

% @author JStiso

% Normalize
if nor == 1
%     A = A/(eigs(A,1)+1) - eye(size(A));
%     max_eig = eigs(A,1);
%     c = 0.01 * max_eig;
%     A = A/(max_eig + c) - eye(size(A));
    max_eig = eigs(A,1);
    A = 0.9*(A/(max_eig)) - eye(size(A));
    disp('Here Mod.')
end

[V, ~] = eig(A);
lambda = eig(A);

% convert to discrete lambda - comment out if system is already discrete
% (mapping to discrete time system)
lambda = exp(lambda);
% get only the X% largest values
[modes, idx] = sort(lambda);
cutoff = round(numel(modes)*thresh);
transient_modes = modes(1:cutoff);
persistent_modes = modes((end-cutoff+1):end);
% get controllability
phi_p = sum((1 - (persistent_modes)).*(V(i,idx((numel(idx)-cutoff+1):end)).^2)');
phi_t = sum((1 - (transient_modes)).*(V(i,idx(1:cutoff)).^2)');

end
