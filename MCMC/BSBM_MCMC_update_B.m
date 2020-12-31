function B = BSBM_MCMC_update_B(m, Gam, l0, u0)

% Simulate B in the MCMC Algorithm
%
% Usage: B = BSBM_MCMC_update_B(m, Gam, l0, u0)
%
% Input: 
% m,Gam         - The values of m,Gam at current iteration
% l0,u0         - Hyperparamters in B's prior
%
% Output:
% B             - The simulated values of B


[K, d] = size(m);

B = zeros(K, K);


% Create the collections of indexes
%--------------------------------------------------------------------------
collect = [];

for i = 1:K
    
        m_tmp = zeros(K, d);
        m_tmp(i, :) = ones(1, d);
        m_tmp = m - m_tmp;
        m_tmp = (m_tmp == 0);
        collect{i} = find(all(m_tmp));
    
end


% Make draws
%--------------------------------------------------------------------------
for i = 1:K

    a1 = 1;
    b1 = 1;
    collect_tmp = collect{i};
    
    H = length(collect_tmp);
    
    for k = 1:H
        for w = 1:H
        
            a1 = a1 + Gam(collect_tmp(k), collect_tmp(w));
            b1 = b1 + 1 - Gam(collect_tmp(k), collect_tmp(w));
            
        end
    end
    
    P1 = betacdf(l0, a1, b1);
    P2 = 1;
    B(i, i) = betainv(P1+(P2-P1)*rand, a1, b1);
    
    for j = 1:K
        if j == i
            continue
        end
    
    collect_j = collect{j};
    W = length(collect_j);
    
    a2 = 1;
    b2 = 1;
    
    for k = 1:H
        for w = 1:W
        
            a2 = a2 + Gam(collect_tmp(k), collect_j(w));
            b2 = b2 + 1 - Gam(collect_tmp(k), collect_j(w));
            
        end
    end
    
        P1 = 0;
        P2 = betacdf(u0, a2, b2);
        B(i, j) = betainv(P1+(P2-P1)*rand, a2, b2);
    
    end
end

end