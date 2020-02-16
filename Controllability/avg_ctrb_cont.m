function v = avg_ctrb_cont(A, T, nor)
% Function for computing average controllability for each node.
% Unsure of numerical stability, don't use for publication results
%
% Inputs:
% A:    N x N continuous-time stable adjacency matrix
% T:    1 x 1 scalar for time horizon
% nor:  normalization constant (SP)
% Outputs:
% v:    N x 1 vector of average controllability values

% Normalize
if nor == 1
%     A = A/(eigs(A,1)+1) - eye(size(A));
%     max_eig = eigs(A,1);
%     c = 0.01 * max_eig;
%     A = A/(max_eig + c) - eye(size(A));
    max_eig = eigs(A,1);
    A = 0.9*(A/(max_eig)) - eye(size(A));
    disp('Here Avg.');
end

% System Parameters
N = size(A,1);

if(T < Inf)
    G = integral(@(t)expm((A+A')*t), 0, T, 'ArrayValued', 1,...
                 'AbsTol', 1e-12, 'RelTol', 1e-12);
    v = diag(G);
elseif(T==Inf)
    sys = ss(A,eye(size(A)),zeros(1,N),0);
    Wc = gram(sys,'c');
    v = diag(Wc);
end

end