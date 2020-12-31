function mu = BSBM_EM_update_mu(mean_T, xi_sq)

% Update mu in M-step of the EM Algorithm
%
% Usage: mu = BSBM_EM_update_mu(mean_T, xi_sq)
%
% Input: 
% mean_T   - The outputs of Kalman Filter and Smoother in the E-step of
%            current iteration
% xi_sq    - Hyperparameter in mu's prior
%
% Output:
% mu       - The updated values of mu


mu = (1/(1+1/xi_sq)).*mean_T(:,1);


end