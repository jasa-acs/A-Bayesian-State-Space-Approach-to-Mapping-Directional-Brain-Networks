function X = BSBM_MCMC_update_X(X_prev, A_gam, C, Q, R, mu, V, Y)

% Simulate X in the MCMC Algorithm
%
% Usage: X = BSBM_MCMC_update_X(X_prev, A_gam, C, Q, R, mu, V, Y)
%
% Input: 
% Y             - The observed data
% X_prev        - The values of X at current iteration before simulating X
% A_gam         - The values of A.*Gam at current iteration
% C,R,mu        - The values of C,R,mu at current iteration
% (Q,V          - Fixed at the identity matrices in the current version of
%                 BSBM)
%
% Output:
% X             - The simulated values of X


[~, T] = size(Y);
d = size(Q, 1);
X = X_prev;

for t = 1:(T+1)
    
% Update x(0)
%--------------------------------------------------------------------------
    if t == 1
        
        U0 = (A_gam' /Q) * A_gam + inv(V);
        U0 = (U0 + U0.')/2;
               
        [dec_U, ~] = cholcov(U0);
        
        V0 = (X(:,2)' / Q) * A_gam + mu' / V;
        mu2 = U0 \ (V0');
        
        tmp = randn(1, d) / (dec_U') + mu2';
        
        X(:, t) = tmp';        
    
        
% Update x(T)
%--------------------------------------------------------------------------
    elseif t == (T+1)

        U0 = inv(Q) + (C'/R) * C;
        V0 = (X(:,t-1)' * A_gam') / Q + (Y(:, t-1)' / R) * C;
                      
        tmp = mvnrnd(U0 \ (V0'), (inv(U0) + (inv(U0)).')/2);
        X(:, t) = tmp';        
        

% Update x(1),x(2),...,x(T-1)
%--------------------------------------------------------------------------
    elseif t == 2 
        
        U_tmp = (A_gam' / Q) * A_gam + inv(Q) + (C'/R) * C;
        inv_U = inv(U_tmp); %inv_U is the same for t = 2,3,...,T, we only calculate it once
        V0 = (X(:,t-1)' * A_gam') / Q + (X(:,t+1)' / Q) * A_gam + (Y(:, t-1)' / R) * C;
        
        tmp = mvnrnd(inv_U * (V0'), (inv_U + inv_U.')/2);
        X(:, t) = tmp';  
        
        
    else

        V0 = (X(:,t-1)' * A_gam') / Q + (X(:,t+1)' / Q) * A_gam + (Y(:, t-1)' / R) * C;
        tmp = mvnrnd(inv_U * (V0'), (inv_U + inv_U.')/2);
        X(:, t) = tmp';  
        
    end
          
end

end