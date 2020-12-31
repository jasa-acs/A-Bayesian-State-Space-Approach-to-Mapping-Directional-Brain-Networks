function [] = BSBM_ContPnl()

% The Control Panel of Our Proposed Bayesian Method with a SBM-motivated
% Prior (BSBM)
% This Control Panel serves as an example of analyzing data using our
% toolbox. It is written as a function which can be easily modified to use
% on computing clusters. Users can use this control panel to reproduce the
% results of our simulation study (Section 3) and modify this file to
% analyze their data without difficulty.


% Add Path
addContainingDirAndSubDir();


%% Data Loading
% Load the data to our workspace. We use the simulated data that used in
% Section 3 as an example. The simulated datasets are provided in the
% subfolder "Simulated_data" in our toolbox:
%
% 'Simulation1.mat'     - The simulated data we used in Section 3.1
% 'Simulation2_1.mat'   - The simulated data we used in Section 3.2, with
%                         2714 time points
% 'Simulation2_2.mat'   - The simulated data we used in Section 3.2, with
%                         1000 time points
%
% The systems and the error structures we used to generate these simulated
% datasets are described in Section 3.1 and Section 3.2 of our paper in
% detail.
% 
% Each data contains two variables:
%
% Y             - The simulated data, which is a real d*T matrix. d is the
%                 number of channels and T is the number of time points
% A_true        - The true network structure. "A_true(i,j) = 1" indicates 
%                 that there is a directional connection from channel j to
%                 channel i


% The data path
filename = "Simulation1.mat";        % We use "Simulation1.mat" as an example here
load(filename); 
% The save path of the outputs of the EM algorithm
filename2 = "Simulation1_EM.mat";    
% The save path of the outputs of the MCMC algorithm
filename3 = "Simulation1_MCMC.mat";  
% The save path of the clustering probabilities and the network edge probabilities
filename4 = "Simulation1_Gam_m.mat"; 


%% Data Standardization

Y = Y';
Y = datastd(Y);
[~,d]=size(Y);
Y = Y';


%% Initialization

rand('seed', 31829);

% Specify the number of iterations for the EM algorithm
mite = 200; 

% Specify the number of iterations for the MCMC algorithm
mcmc_mite = 10000;

% Specify hyperparameters for B's prior
l0 = 0.9;   
u0 = 0.1;

% Specify the hyperparameter for R's prior
e0 = 0.000001;

% Specify the hyperparameter for A's prior
xi_sq = 0.5;

% Specify the hyperparameter for C's prior and mu's prior
xi_sq2 = 100;

b_sq = 1; % Fixed at 1 in current version of BSBM

K = d;

C0 = eye(d);

mu0 = Y(:,1);

R0 = eye(d);

m0 = eye(d);

B0 = eye(K) * l0 + (ones(K, K) - eye(K)) * u0;

Gam0 = binornd(1, 0.5, d, d);

A0 = eye(d);


%% Perform EM Algorithm

tic

[Gam_all, A_all, C_all, R_all, mu_all, m_all, B_all] = ...
    BSBM_EM(Y, mite, l0, u0, e0, xi_sq, xi_sq2, A0, C0, mu0, R0, m0, B0, Gam0, b_sq);

toc

save(filename2);


%% Set the Outputs of EM as the Starting Values for MCMC

% % Specify the number of iterations for the MCMC algorithm
% mcmc_mite = 10000;

m0 = m_all(:,:, mite+1);

tmp_k = zeros(d,1);
for i = 1:d
    tmp_k(i) = find(m0(:,i));
end
tmp_k2 = unique(tmp_k, "stable");
tmp_l = length(tmp_k2);
if tmp_l > 1
    mcmc_K = tmp_l;
else
    mcmc_K = 4;
end

m0 = m0(1:mcmc_K, :);

C0 = C_all(:,:,mite+1);

R0 = R_all(:,:,mite+1);

mu0 = mu_all(:,mite+1);

B0 = B_all(:,:,mite+1);
B0 = B0(1:mcmc_K, 1:mcmc_K);

A0 = A_all(:,:,mite+1);

Gam0 = Gam_all(:,:,mite+1);

Q0 = eye(d);
V0 = eye(d);

tmp_m = [m0 10^-10*ones(mcmc_K,1)];
P0 = sum(tmp_m, 2)/sum(tmp_m, [1 2]);
P0 = P0.*ones(mcmc_K,d);

X0 = horzcat(mu0, Y);


tic

[A_all, C_all, X, V, Q, R_all, mu_all, P_all, m_all, B_all, Gam_all] ...
    = BSBM_MCMC(Y, mcmc_mite, l0, u0, e0, xi_sq, xi_sq2, X0, A0, C0, mu0, Q0, ...
    R0, V0, P0, m0, B0, Gam0);

toc

save(filename3);


%% Calculate the Clustering Probabilities and the Network Edge Probabilities
% We use the clustering probabilities and the network edge probabilities to
% map the brain network:
%
% m_result     - The clustering probabilities
% Gam_result   - The network edge probabilities
% 
% The analyses based on the clustering probabilities and the network edge
% probabilities are described in Section 2.4, Section 3, and Section 4 of
% our paper in detail.


m_result = zeros(d,d);
for nr = 1:d
    for nc = nr:d
        z = 0;
        for s = 2002:10001
            if m_all(:,nr,s) == m_all(:,nc,s)
                z = z+1;
            end
        end
        m_result(nr,nc) = z/8000;
    end
end
for nr2 = 1:d
    for nc2 = 1:nr2
        m_result(nr2, nc2) = m_result(nc2,nr2);
    end
end

Gam_result = mean(Gam_all(:,:,2002:10001), 3);

save(filename4, 'Gam_result', 'm_result', 'A_true');

%save(filename4, 'Gam_result', 'm_result'); % For analyzing real data

end