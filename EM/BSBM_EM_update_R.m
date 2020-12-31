function R = BSBM_EM_update_R(mean_T, var_T, Y, C, e0)

% Update R in M-step of the EM Algorithm
%
% Usage: R = BSBM_EM_update_R(mean_T, var_T, Y, C, e0)
%
% Input: 
% Y              - The observed data
% mean_T,var_T   - The outputs of Kalman Filter and Smoother in the E-step
%                  of current iteration
% C              - The values of C at current iteration before updating R
% e0             - Hyperparameter in R's prior
%
% Output:
% R              - The updated values of R


[d, T] = size(Y);
R = eye(d);

Phi = zeros(d,d);
Psi = zeros(d,T);

for t = 1:T
    Phi = Phi + C*var_T(:,:,(t+1))*C;
    Psi(:,t) = Y(:,t) - C*mean_T(:,(t+1));
end

for i = 1:d
    R(i,i) = (e0 + Phi(i,i)/2 + Psi(i,:)*Psi(i,:).'/2)/(T/2 + 1 + e0);
end

end