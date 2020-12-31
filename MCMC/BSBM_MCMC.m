function [A_all, C_all, X, V, Q, R_all, mu_all, P_all, m_all, B_all, Gam_all] ...
    = BSBM_MCMC(Y, mite, l0, u0, e0, xi_sq, xi_sq2, X0, A0, C0, mu0, Q0, R0, V0, P0, m0, B0, Gam0)

%% Description
% Partial Collapsed Gibbs Sampler to provide simulations from the posterior
% distribution
%
% Usage:
% [A_all, C_all, X, V, Q, R_all, mu_all, P_all, m_all, B_all, Gam_all] =
% BSBM_MCMC(Y, mite, l0, u0, e0, xi_sq, xi_sq2, X0, A0, C0, mu0, Q0, R0, V0, P0,
% m0, B0, Gam0)
%
% Input:
% Y                            - The observed data 
% mite                         - The max number of iterations
% l0,u0                        - Hyperparamters in B's prior
% e0                           - Hyperparameter in R's prior
% xi_sq                        - Hyperparameter in A's prior
% xi_sq2                       - Hyperparameter in C's prior and mu's prior
% X0                           - The initial values of the underlying
%                                neuronal states
% A0,C0,mu0,R0,P0,m0,B0,Gam0   - The initial values of A,C,mu,R,P,m,B,Gam
% (Q0,V0                       - Fixed at the identity matrices in current
%                                version of BSBM)
% 
% Output:
% Gam_all,A_all,C_all,R_all,mu_all,P_all,m_all,B_all record the simulated
% values of Gam,A,C,R,mu,P,m,B at each iteration. X record the simulated
% values of X at the last iteration. 
% (V,Q are fixed at the identity matrices in the current version of BSBM)


%% Initialization

[K, d] = size(P0);

ite = 1;

A_all = zeros(d, d, mite+1);
A_all(:,:,1) = A0;

C_all = zeros(d, d, mite+1);
C_all(:,:,1) = C0;

mu_all = zeros(d, mite+1);
mu_all(:,1) = mu0;

Q = Q0;
V = V0;

R_all = zeros(d, d, mite+1);
R_all(:,:,1) = R0;

X = X0;

P_all = zeros(K, d, mite+1);
P_all(:,:,1) = P0;

m_all = zeros(K, d, mite+1);
m_all(:,:,1) = m0;

B_all = zeros(K, K, mite+1);
B_all(:,:,1) = B0;

Gam_all = zeros(d, d, mite+1);
Gam_all(:,:,1) = Gam0;


%% Iterations

while ite <= mite
    
disp(ite);

% Update Gamma

Gam_all(:,:,ite+1) = BSBM_MCMC_update_Gam(Gam_all(:,:,ite), xi_sq, Q, X, ...
    B_all(:,:,ite), m_all(:,:,ite));

% Update A    
       
A_all(:,:,ite+1) = BSBM_MCMC_update_A(X, Q, Gam_all(:,:,ite+1), xi_sq);
    
% Update X

A_gam = A_all(:,:,ite+1) .* Gam_all(:,:,ite+1);

X = BSBM_MCMC_update_X(X, A_gam, C_all(:,:,ite), Q, R_all(:,:,ite), mu_all(:,ite), V, Y);
   
% Update C

C_all(:,:,ite+1) = BSBM_MCMC_update_C(X, Y, R_all(:,:,ite), xi_sq2);
 
% Update Q (we fix Q at eye(d))

Q = BSBM_MCMC_update_Q(X, A_all(:,:,ite+1), Gam_all(:,:,ite+1));   
       
% Update R 

R_all(:,:,ite+1) = BSBM_MCMC_update_R(X, Y, C_all(:,:,ite+1), e0);

% Update mu

mu_all(:,ite+1) = BSBM_MCMC_update_mu(X, V, xi_sq2);

% Update P

P_all(:,:,ite+1) = BSBM_MCMC_update_P(m_all(:,:,ite));

% Update m

m_all(:,:,ite+1) = BSBM_MCMC_update_m(m_all(:,:,ite), B_all(:,:,ite), Gam_all(:,:,ite+1), P_all(:,:,ite+1));

% Update B

B_all(:,:,ite+1) = BSBM_MCMC_update_B(m_all(:,:,ite+1), Gam_all(:,:,ite+1), l0, u0);

% Update iteration number

ite = ite + 1;

end
end
