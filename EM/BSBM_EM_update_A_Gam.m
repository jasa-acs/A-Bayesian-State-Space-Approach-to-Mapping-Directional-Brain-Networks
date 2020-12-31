function [A_i, Gam_i] = BSBM_EM_update_A_Gam(mean_T, var_T, cov_T, A_prev, ...
    Gam_prev, m, B, xi_sq, b_sq, i)

% Jointly Update A and Gam in M-step of the EM Algorithm
% This function jointly update a specified row of A and Gam.
%
% Usage: [A_i, Gam_i] = BSBM_EM_update_A_Gam(mean_T, var_T, cov_T, A_prev,
% Gam_prev, m, B, xi_sq, b_sq, i)
%
% Input: 
% mean_T,var_T,cov_T    - The outputs of Kalman Filter and Smoother in the 
%                         E-step of current iteration
% i                     - The specified row number i
% A_prev, Gam_prev      - The values of A,Gam at current iteration before
%                         updating the ith row of A,Gam
% xi_sq                 - Hyperparameter in A's prior
% (b_sq                 - Fixed at 1 in current version of BSBM)
%
% Output:
% A_i, Gam_i            - The updated values of the ith row of A,Gam


[d, ~] = size(mean_T);
[~,~,T] = size(cov_T);

A = A_prev;
Gam = Gam_prev;

Upsilon = zeros(d,d);
Delta = zeros(d,d);
for t = 1:T
    Upsilon = Upsilon + cov_T(:,:,t);
    Delta = Delta + var_T(:,:,t);
end


for j = 1:d
        max1 = log(1-m(:,i).'*B*m(:,j)); 
        L1 = Delta(j,j) + 1/(b_sq*xi_sq);
        for t = 1:T
            L1 = L1 + mean_T(j,t)^2;
        end
        L2 = Upsilon(i,j);
        for t = 1:T
            L2 = L2 + mean_T(i, t+1)*mean_T(j,t);
        end
        for k = 1:d
            if k == j
                continue
            end
            if Gam(i,k) == 0
                continue
            end
            L2 = L2 - Delta(k,j)*A(i,k);
            for t = 1:T
                L2 = L2 - mean_T(k,t)*mean_T(j,t)*A(i,k);
            end
        end
        max2 = log(m(:,i).'*B*m(:,j)) + (L2^2)/(2*L1);
        
        if max2 > max1
            Gam(i,j) = 1;
            A(i,j) = L2/L1;
        else
            Gam(i,j) = 0;
            A(i,j) = 0;
        end
end

A_i = A(i,:);
Gam_i = Gam(i,:);

end