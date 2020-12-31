function m = BSBM_MCMC_update_m(m_prev, B, Gam, P)

% Simulate m in the MCMC Algorithm
%
% Usage: m = BSBM_MCMC_update_m(m_prev, B, Gam, P)
%
% Input: 
% m_prev        - The values of m at current iteration before simulating m
% B,Gam,P       - The values of B,Gam,P at current iteration
%
% Output:
% m             - The simulated values of m

     
[K, d] = size(m_prev);
m = m_prev;

for i = 1:d
    J = ones(1, K);
    for l = 1:K
        m_l = zeros(K, 1);
        m_l(l) = 1;
        m_tmp = m;
        m_tmp(:, i) = m_l;
        for j = 1:d 
            if j ~= i
            J(l) = J(l) * ((B(l, :) * m_tmp(:, j))^Gam(i, j)) ...
                * ((1 - B(l, :) * m_tmp(:, j))^(1 - Gam(i, j)))...
                * (((m_tmp(:,j))'*B(:, l))^Gam(j, i)) ...
                * ((1-(m_tmp(:,j))'*B(:, l))^(1-Gam(j, i)));
            end
        end
        J(l) = J(l) * ((B(l, l))^Gam(i, i)) * ((1 - B(l, l))^(1 - Gam(i, i))) * P(l, i);            
    end
    if ~any(J) == 1
        J1 = ones(1, K) / K;
    else
        J1 = J / sum(J);
    end
    
    m_tmp = mnrnd(1, J1, 1);
    m(:, i) = m_tmp';
    
end
end