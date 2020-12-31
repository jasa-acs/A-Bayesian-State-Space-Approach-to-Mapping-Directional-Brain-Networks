function m_i = BSBM_EM_update_m(m_prev, Gam, B, i)

% Update m in M-step of the EM Algorithm
% This function update a specified column of m.
%
% Usage: m_i = BSBM_EM_update_m(m_prev, Gam, B, i)
%
% Input: 
% i              - The specified column number i
% m_prev,Gam,B   - The values of m,Gam,B at current iteration before
%                  updating the ith column of m
%
% Output:
% m_i            - The updated values of the ith column of m


[K, d] = size(m_prev);
m = m_prev;
alpha = sum(m_prev, 2) + ones(K,1);
I = eye(K);

    tmp1 = zeros(K,1);
    for k = 1:K
        tmp_m = I(:,k);
        tmp2 = 0;
        for j = 1:d
            if j == i
                continue
            end
            
            if Gam(i,j) == 1
                tmp2 = tmp2 + log(tmp_m.'*B*m(:,j));
            else 
                tmp2 = tmp2 + log(1-tmp_m.'*B*m(:,j));
            end
            if Gam(j,i) == 1
                tmp2 = tmp2 + log(m(:,j).'*B*tmp_m);
            else 
                tmp2 = tmp2 + log(1-m(:,j).'*B*tmp_m);
            end

        end
        if Gam(i,i) == 1
            tmp2 = tmp2 + log(tmp_m.'*B*tmp_m);
        else
            tmp2 = tmp2 + log(1 - tmp_m.'*B*tmp_m);
        end
        
        %Integrate P out(expectation)
        tmp4 = 0;
        for l = 1:K
            tmp4 = tmp4 + tmp_m(l)*(psi(alpha(l)) - psi(sum(alpha)));
        end
        
        tmp2 = tmp2 + tmp4;
        
        tmp1(k,1) = tmp2;
       
    end
    [~, idx] = max(tmp1);
    m(:,i) = I(:,idx);
    

m_i = m(:,i);



end