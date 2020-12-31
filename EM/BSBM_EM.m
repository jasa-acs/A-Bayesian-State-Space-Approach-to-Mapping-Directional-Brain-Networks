function [Gam_all, A_all, C_all, R_all, mu_all, m_all, B_all] ...
   = BSBM_EM(Y, mite, l0, u0, e0, xi_sq, xi_sq2, A0, C0, mu0, R0, m0, B0, Gam0, b_sq)

%% Description
% EM Algorithm for Setting Starting Values and Hyperparameter of the
% Partially Collapsed Gibbs Sampler
%
% Usage:
% [Gam_all, A_all, C_all, R_all, mu_all, m_all, B_all] = BSBM_EM(Y, mite, l0,
% u0, e0, xi_sq, xi_sq2, A0, C0, mu0, R0, m0, B0, Gam0, b_sq)
% 
% Input:
% Y                         - The observed data 
% mite                      - The max number of iterations
% l0,u0                     - Hyperparamters in B's prior
% e0                        - Hyperparameter in R's prior
% xi_sq                     - Hyperparameter in A's prior
% xi_sq2                    - Hyperparameter in C's prior and mu's prior
% A0,C0,mu0,R0,m0,B0,Gam0   - The initial values of A,C,mu,R,m,B,Gam
% (b_sq                     - Fixed at 1 in current version of BSBM)
%
% Output:
% Gam_all,A_all,C_all,R_all,mu_all,m_all,B_all record the updated values of
% Gam,A,C,R,mu,m,B at each iteration. 


%% Initialization
[K, d] = size(m0);

ite = 1;

A_all = zeros(d, d, mite+1);
A_all(:,:,1) = A0;
C_all = zeros(d, d, mite+1);
C_all(:,:,1) = C0;
mu_all = zeros(d, mite+1);
mu_all(:,1) = mu0;
R_all = zeros(d, d, mite+1);
R_all(:,:,1) = R0;
m_all = zeros(K, d, mite+1);
m_all(:,:,1) = m0;
B_all = zeros(K, K, mite+1);
B_all(:,:,1) = B0;
Gam_all = zeros(d, d, mite+1);
Gam_all(:,:,1) = Gam0;


%% Iterations

while ite <= mite
    
disp(ite);

% E-step:
%--------------------------------------------------------------------------
[mean_T, var_T, cov_T] = BSBM_Kalman(mu_all(:,ite), Gam_all(:,:,ite),...
    A_all(:,:,ite), C_all(:,:,ite), R_all(:,:,ite), Y);


% M-step:
%--------------------------------------------------------------------------
% Update C
C_all(:,:,ite+1) = BSBM_EM_update_C(mean_T, var_T, Y, R_all(:,:,ite), xi_sq2); 

% Update R
R_all(:,:,ite+1) = BSBM_EM_update_R(mean_T, var_T, Y, C_all(:,:,ite+1), e0);

% Update mu
mu_all(:,ite+1) = BSBM_EM_update_mu(mean_T, xi_sq2);

% Update A, Gam, m
tmp_A = A_all(:,:,ite);
tmp_Gam = Gam_all(:,:,ite);
tmp_m = m_all(:,:,ite);

for node_i = 1:d

tmp_m(:,node_i) = BSBM_EM_update_m(tmp_m, tmp_Gam, B_all(:,:,ite), node_i);   
[tmp_A(node_i,:), tmp_Gam(node_i,:)] = BSBM_EM_update_A_Gam(mean_T, var_T, ...
    cov_T, tmp_A, tmp_Gam, tmp_m, B_all(:,:,ite), xi_sq, b_sq, node_i);

end
A_all(:,:,ite+1) = tmp_A;
Gam_all(:,:,ite+1) = tmp_Gam;
m_all(:,:,ite+1) = tmp_m;

% Update B
B_all(:,:,ite+1)= BSBM_EM_update_B(m_all(:,:,ite+1), Gam_all(:,:,ite+1), B_all(:,:,ite), l0, u0);


% Sort the rows of m,B
%--------------------------------------------------------------------------
% Sort the rows of m
m = m_all(:,:,ite+1);
tmp_k = zeros(d,1);
for i = 1:d
    tmp_k(i) = find(m(:,i));
end
tmp_k2 = unique(tmp_k, "stable");
m2 = zeros(K,d);
for i = 1:d
    m2(tmp_k2 == tmp_k(i),i) = 1;
end

% Sort the rows of B
B = B_all(:,:,ite+1);
tmp_diff = setdiff(1:K, tmp_k2);
tmp_order = [tmp_k2;tmp_diff.'];
B2 = B(tmp_order, tmp_order);

m_all(:,:,ite+1) = m2;
B_all(:,:,ite+1) = B2;


% Update the iteration number
%--------------------------------------------------------------------------
ite = ite + 1;

end
end