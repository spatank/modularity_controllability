function Q = modularity_index(A,Ci)
% modularity_index Computes the modularity quality index Q given 
%                  network partitions defined by Ci
%   A: network adjacency matrix
%   Ci: community assignments

N = length(A); 
K = sum(A);
m = sum(K);
gamma = 1; % classic modularity
B = A - gamma*(K.'*K)/m;
s = Ci(:,ones(1,N)); % compute modularity
Q = ~(s-s.').*B/m;
Q = sum(Q(:));
end

