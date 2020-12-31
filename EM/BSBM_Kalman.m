function [mean_T, var_T, cov_T] = BSBM_Kalman(mu, Gam, A, C, R, Y)

% Kalman Filter and Kalman Smoother in E-step of the EM algorithm
%
% Usage: [mean_T, var_T, cov_T] = BSBM_Kalman(mu, Gam, A, C, R, Y)
% 
% Input:
% mu,Gam,A,C,R,Y       - The values of mu,Gam,A,C,R,Y at current iteration
%                        (the updated values from the last iteration)
%
% Output:
% mean_T,var_T,cov_T   - Components of the log posterior probability after
%                        integrating the underlying neuronal states out


[~, T] = size(Y);
d = size(R, 1);

A_gam = A.*Gam;

%% Kalman Filter
mean_t = zeros(d, (T+1));
var_t = zeros(d, d, (T+1));

mean_t(:, 1) = mu;
var_t(:, :, 1) = eye(d);

mean_t_1 = zeros(d, (T+1));
var_t_1 = zeros(d, d, (T+1));
K_t = zeros(d,d, (T+1));

for t = 2:(T+1)
    mean_t_1(:,t) = A_gam*mean_t(:, (t-1));
    var_t_1(:,:,t) = A_gam*var_t(:,:, (t-1))*A_gam.' + eye(d);
    var_t_1(:,:,t) = (var_t_1(:,:,t) + var_t_1(:,:,t).')/2;
    inv_tmp = C*var_t_1(:,:,t)*C + R;
    inv_tmp = (inv_tmp + inv_tmp.')/2;
    
    K_t(:,:,t) = (var_t_1(:,:,t)*C)/inv_tmp;
    mean_t(:, t) = mean_t_1(:,t) + K_t(:,:,t)*(Y(:,(t-1)) - C*mean_t_1(:,t));
    var_t(:, :, t) = (eye(d) - K_t(:,:,t)*C)*var_t_1(:,:,t);
    var_t(:, :, t) = (var_t(:, :, t) + var_t(:, :, t).')/2;
end


%% Kalman Smoother
mean_T = zeros(d, (T+1));
var_T = zeros(d,d, (T+1));

mean_T(:,(T+1)) = mean_t(:,(T+1));
var_T(:,:,(T+1)) = var_t(:,:,(T+1));

J_t = zeros(d,d,(T+1));

for i  = 1:T
    t = (T+1) - i;
    J_t(:,:,t) = (var_t(:,:,t)*A_gam.')/var_t_1(:,:,(t+1));
    mean_T(:,t) = mean_t(:,t) + J_t(:,:,t)*(mean_T(:,(t+1)) - mean_t_1(:,(t+1)));
    var_T(:,:,t) = var_t(:,:,t) + J_t(:,:,t)*(var_T(:,:,(t+1)) - var_t_1(:,:,(t+1)))*J_t(:,:,t).';
    var_T(:,:,t) = (var_T(:,:,t) + var_T(:,:,t).')/2;
end


%% Lag-One Covariance Smoother
cov_T = zeros(d,d,T);
cov_T(:,:,T) = (eye(d) - K_t(:,:,(T+1))*C)*A_gam*var_t(:,:,T);

for i = 1:(T-1)
    t = T-i;
    cov_T(:,:,t) = var_t(:,:,(t+1))*J_t(:,:,t).' +...
        J_t(:,:,(t+1))*(cov_T(:,:,(t+1))-A_gam*var_t(:,:,(t+1)))*J_t(:,:,t).';
end


end