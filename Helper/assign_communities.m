function [labels] = assign_communities(Model)
% assign_communities takes in a fitted WSBM and assigns community labels to
% the nodes in the network
%   Model: fitted Weighted Stochastic Block Model

mus = Model.Para.mu;

labels = zeros(length(mus),1); 

for i = 1:length(labels)
    [~, labels(i)] = max(mus(:,i)) ; 
end
    
end
