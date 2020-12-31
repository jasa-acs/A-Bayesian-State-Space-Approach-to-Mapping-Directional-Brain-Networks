function C = BSBM_MCMC_update_C(X, Y, R, xi_sq)

% Simulate C in the MCMC Algorithm
%
% Usage: C = BSBM_MCMC_update_C(X, Y, R, xi_sq)
%
% Input: 
% Y             - The observed data
% X             - The values of X at current iteration
% R             - The values of R at current iteration
% xi_sq         - Hyperparameter in C's prior
%
% Output:
% C             - The simulated values of C

 
[d, T] = size(Y);
C = eye(d);


for i = 1:d
    
U_i = 0;
V_i = 0;

for t = 1:T
    
    U_i = U_i + (X(i, t+1))^2;
    V_i = V_i + X(i, t+1) * Y(i, t);
  
end
    
U_i = U_i / R(i, i) + 1/xi_sq;
V_i = V_i / R(i, i);

sigma = 1 / U_i;
mu = V_i / U_i;
    
    
C(i, i) = randn * sqrt(sigma) + mu;


end

end
