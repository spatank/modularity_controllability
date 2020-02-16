function v = avg_ctrb_disc(A, T, nor)
% Function for computing discrete average controllability for each node.
% Unsure of numerical stability, don't use for publication results
% Infinite horizon code from Shi Gu et al. 2015 doi: 10.1038/ncomms9414
%
% Inputs:
% A:    N x N continuous-time stable adjacency matrix
% T:    1 x 1 scalar integer for time horizon (>0)
% nor:  normalization constant (default: 1+eig_largest) (SP)
% Outputs:
% v:    N x 1 vector of average controllability values

% Normalize
% if nor == 1
%     A = A./(1+eigs(A,1));
% end
if nargin ~= 3
    nor = (1+eigs(A,1));
    A = A./nor;
else
    A = A./nor;
end

% System Parameters
N = size(A,1);

if(T < Inf)
    G = zeros(N);
    for i = 0:T
        G = G + A^(2*i);
    end
    v = diag(G);
elseif(T==Inf)
    [U, T] = schur(A,'real'); % Schur stability
    midMat = (U.^2)';
    v = diag(T);
    P = repmat(diag(1 - v*v'),1,size(A,1));
    v = sum(midMat./P)';
end

end