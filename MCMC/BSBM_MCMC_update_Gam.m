function Gam = BSBM_MCMC_update_Gam(Gam_prev, xi_sq, Q, X, B, m)

% Simulate Gam in the MCMC Algorithm
%
% Usage: Gam = BSBM_MCMC_update_Gamma(Gam_prev, xi_sq, Q, X, B, m)
%
% Input: 
% Gam_prev      - The values of Gam at current iteration before simulating
%                 Gam
% X             - The values of X at current iteration
% B,m           - The values of B,m at current iteration
% (Q            - Fixed at the identity matrix in the current version of
%                 BSBM)
% xi_sq         - Hyperparameter in A's prior
%
% Output:
% Gam           - The simulated values of Gam

    
[d, T] = size(X);
T = T-1;
Gam = Gam_prev;

for i = 1:d
    
    for j = 1:d
        
        Gamij1 = Gam(i, :);
        Gamij0 = Gam(i, :);        
        Gamij1(j) = 1;
        Gamij0(j) = 0;        
        
        X1_tmp = Gamij1' .* X(:,1:T);
        X0_tmp = Gamij0' .* X(:,1:T);
        U1 = X1_tmp * X1_tmp';
        U0 = X0_tmp * X0_tmp';
        V1_tmp = sum(X(i, 2:(T+1))' .* X1_tmp', 1);
        V0_tmp = sum(X(i, 2:(T+1))' .* X0_tmp', 1);
        V1_tmp = V1_tmp / Q(i,i);
        V0_tmp = V0_tmp / Q(i,i);
        U1 = U1 / Q(i,i) + eye(d)/xi_sq;
        U0 = U0 / Q(i,i) + eye(d)/xi_sq;
        
        invU0 = inv(U0);
        invU1 = inv(U1);
        
        logw1 = 0.5*log(det(invU1)) + log((m(:,i)' * B * m(:,j))) + 0.5 * V1_tmp * (invU1 * V1_tmp'); 
        logw0 = 0.5*log(det(invU0)) + log(1 - m(:,i)' * B * m(:,j)) + 0.5 * V0_tmp * (invU0 * V0_tmp'); 
                
        if m(:,i)' * B * m(:,j) == 1
            
           prob = 1;           
                 
        elseif m(:,i)' * B * m(:,j) == 0
            
           prob = 0;
                 
        elseif logw0 == -Inf && logw1 == -Inf
            
           prob = 0.5; 
           
        elseif logw0 == Inf && logw1 == Inf
            
           prob = 0.5;    
           
        else    
            
            prob = 1 / (exp(logw0 - logw1) + 1);
                      
        end
        
        Gam(i, j) = binornd(1, prob);
        
                 
    end

    
end


end
