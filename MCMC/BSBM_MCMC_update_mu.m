function mu = BSBM_MCMC_update_mu(X, V, xi_sq)

% Simulate mu in the MCMC Algorithm
%
% Usage: mu = BSBM_MCMC_update_mu(X, V, xi_sq)
%
% Input: 
% X             - The values of X at current iteration
% xi_sq         - Hyperparameter in mu's prior
%
% Output:
% mu            - The simulated values of mu


[d, ~] = size(X);

sigma = inv(inv(V) + (1/xi_sq)*eye(d));

tmp = sigma * (inv(V) * X(:,1));

mu = mvnrnd(tmp', sigma);    

mu = mu';

end