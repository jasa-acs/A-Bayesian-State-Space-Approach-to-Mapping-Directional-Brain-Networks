function P = BSBM_MCMC_update_P(m)

% Simulate P in the MCMC Algorithm
%
% Usage: P = BSBM_MCMC_update_P(m)
%
% Input: 
% m             - The values of m at current iteration
%
% Output:
% P             - The simulated values of P


[K, d] = size(m);

m_sum = sum(m, 2);

P = zeros(K, d);

a = m_sum + ones(K, 1);
a = a';

P_tmp = gamrnd(repmat(a,1,1),1,1,length(a));

P_tmp = P_tmp ./ repmat(sum(P_tmp,2),1,length(a));

P_tmp = P_tmp';

P(:, 1) = P_tmp;

for i = 2:d
    
P(:,i) = P(:, 1);
    
end

end


