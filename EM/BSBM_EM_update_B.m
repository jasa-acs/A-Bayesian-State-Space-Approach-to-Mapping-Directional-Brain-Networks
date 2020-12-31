function B = BSBM_EM_update_B(m, Gam, B_prev, l0, u0)

% Update B in M-step of the EM Algorithm
%
% Usage: B = BSBM_EM_update_B(m, Gam, B_prev, l0, u0)
%
% Input: 
% m,Gam,B_prev   - The values of m,Gam,B at current iteration before
%                  updating B
% l0,u0          - Hyperparamters in B's prior
%
% Output:
% B              - The updated values of B


[K, d] = size(m);
B = B_prev;

% Update diagonal entries of B
for k = 1:K
    tmp1 = 0;
    tmp2 = 0;
    for i = 1:d
        for j = 1:d
            if (m(k,i)==1) && (m(k,j)==1)
                tmp1 = tmp1 + 1;
                tmp2 = tmp2 + Gam(i,j);
            end
        end
    end
    
    if tmp1 == 0
        continue
    else
        B(k,k) = max([tmp2/tmp1, l0]);
    end
    
    if B(k,k) == 1
        B(k,k) = 1 - 10^(-10);
    end
    
    if B(k,k) == 0
        B(k,k) = 10^(-10);
    end
end

% Update off-diagonal entries of B
for g = 1:K
    for h = 1:K
        if h == g
            continue
        end
        tmp1 = 0;
        tmp2 = 0;
        for i = 1:d
            for j = 1:d
                if (m(g,i)==1) && (m(h,j)==1)
                    tmp1 = tmp1 + 1;
                    tmp2 = tmp2 + Gam(i,j);
                end
            end
        end
        if tmp1 == 0
            continue
        else
            B(g,h) = min([tmp2/tmp1, u0]);
        end
        if B(g,h) == 1
            B(g,h) = 1 - 10^(-10);
        end
        
        if B(g,h) == 0
            B(g,h) = 10^(-10);
        end
    end
end


end