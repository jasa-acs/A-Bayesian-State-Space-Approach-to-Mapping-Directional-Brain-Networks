function R = BSBM_MCMC_update_R(X, Y, C, e0)

% Simulate R in the MCMC Algorithm
%
% Usage: R = BSBM_MCMC_update_R(X, Y, C, e0)
%
% Input: 
% Y             - The observed data
% X             - The values of X at current iteration
% C             - The values of C at current iteration
% e0            - Hyperparameter in R's prior
%
% Output:
% R             - The simulated values of R


[d, T] = size(Y);
R = eye(d);

for i = 1:d
   
    a = T/2 + e0;
    b = 0;
    
    for t = 1:T
        
        b_tmp = Y(:, t) - C * X(:, t+1);
        b_tmp = b_tmp(i);
        b_tmp = b_tmp^2;
        b = b + b_tmp;  
        
    end
    
    b = b/2 + e0;
        
    R(i, i) = 1 / gamrnd(a, 1/b);

    
end
end