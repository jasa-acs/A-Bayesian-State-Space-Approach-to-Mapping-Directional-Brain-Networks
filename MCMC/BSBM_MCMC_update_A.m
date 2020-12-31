function A = BSBM_MCMC_update_A(X, Q, Gam, xi_sq)

% Simulate A in the MCMC Algorithm
%
% Usage: A = BSBM_MCMC_update_A(X, Q, Gam, xi_sq)
%
% Input: 
% X             - The values of X at current iteration
% Gam           - The values of Gam at current iteration
% (Q            - Fixed at the identity matrix in the current version of
%                 BSBM)
% xi_sq         - Hyperparameter in A's prior
%
% Output:
% A             - The simulated values of A

 
[d, T] = size(X);

T = T-1;

A = zeros(d, d);

for i = 1:d
    
    X_tmp = Gam(i, :)' .* X(:,1:T);
    U_i = X_tmp * X_tmp';
    V_i = sum(X(i, 2:(T+1))' .* X_tmp', 1);
    mu = (U_i / Q(i, i) + (1/xi_sq) * eye(d)) \ (V_i / Q(i, i))';
    U_i = U_i / Q(i, i) + (1/xi_sq) * eye(d);
    
    
    [dec_U, ~] = cholcov(U_i);
        
    A(i, :) = randn(1, d) / (dec_U') + mu';
        

end

end


