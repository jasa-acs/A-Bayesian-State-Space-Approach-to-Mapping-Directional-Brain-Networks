function C = BSBM_EM_update_C(mean_T, var_T, Y, R, xi_sq)

% Update C in M-step of the EM Algorithm
%
% Usage: C = BSBM_EM_update_C(mean_T, var_T, Y, R, xi_sq)
%
% Input: 
% Y              - The observed data
% mean_T,var_T   - The outputs of Kalman Filter and Smoother in the E-step
%                  of current iteration
% R              - The values of R at current iteration before updating C
% xi_sq          - Hyperparameter in C's prior
%
% Output:
% C              - The updated values of C


[d,T] = size(Y);
C = eye(d);

Lambda = zeros(d,d);
for t = 1:T
    Lambda = Lambda + var_T(:,:,(t+1));
end

for i = 1:d
    xi = mean_T(i,2:(T+1)).';
    yi = Y(i,:).';
    C(i,i) = (xi.'*yi)/(Lambda(i,i) + R(i,i)/xi_sq + xi.'*xi);
end

end