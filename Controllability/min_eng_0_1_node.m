function v = min_eng_0_1_node(A, T, sys)
% Function for computing minimum control energy for state transition to 
% unit norm random final states.
% Unsure of numerical stability, don't use for publication results
%
% Inputs:
% A:    N x N continuous-time stable adjacency matrix
% T:    1 x 1 scalar for time horizon
% nor:  normalization constant (SP)
% Outputs:
% v:    N x 1 vector of energy values

% System Parameters
N = size(A,1);

if sys == 'cont'
    if(T < Inf)
        G = integral(@(t)expm(A*t)*expm(A'*t), 0, T, 'ArrayValued', 1,...
                     'AbsTol', 1e-12, 'RelTol', 1e-12);
        G_inv = inv(G);
        v = diag(G_inv);
    elseif(T==Inf)
        sys = ss(A,eye(size(A)),zeros(1,N),0);
        Wc = gram(sys,'c');
        Wc_inv = inv(Wc);
        v = diag(Wc_inv);
    end
end

if sys == 'disc'
    if(T < Inf)
        G = zeros(N);
        for i = 0:T
            G = G + A^(2*i);
        end
        v = diag(inv(G));
    elseif(T==Inf)
        [U, T] = schur(A,'real'); % Schur stability
        midMat = (U.^2)';
        v = diag(T);
        P = repmat(diag(1 - v*v'),1,size(A,1));
        v = 1./(sum(midMat./P)');
    end
end

end