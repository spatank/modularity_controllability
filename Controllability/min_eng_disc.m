function [ x, u, n_err ] = min_eng_disc(A, T, B, x0, xf, nor)
% Computes minimum control energy for state transition.
% A: System adjacency matrix:       N x N
% B: Control input matrix:          N x k
% x0: Initial state:                N x 1
% xf: Final state:                  N x 1
% T: Control horizon                1 x 1
% 
% Outputs
% x: State Trajectory               T+1 x N
% u: Control Input                  T x k

% Normalize
if nor == 1
    A = A/(eigs(A,1)+1);
end 

% System Size
n = size(A,1);
k = size(B,2);

% Compute Controllability Matrix
C = zeros(n,k*T);
for i = 0:T-1
    C(:,[1:k]+(k*i)) = A^i*B;
end

% Compute Control Input
u = pinv(C)*(xf - A^T*x0);
u = reshape(u, [k,T]);
u = fliplr(u);

% Simulate Trajectory
x = zeros(n,T+1);
x(:,1) = x0;
for i = 2:T+1
    x(:,i) = A*x(:,i-1) + B*u(:,i-1);
end

% Error
n_err = norm(xf - x(:,end));
disp(n_err);

% transpose to be similar to opt_eng_cont
u = u';
x = x';
end