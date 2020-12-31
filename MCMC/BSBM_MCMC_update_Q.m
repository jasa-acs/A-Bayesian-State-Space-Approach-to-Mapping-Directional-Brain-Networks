function Q = BSBM_MCMC_update_Q(X, A, Gam)

% We fix Q at the identity matrix in the current version of BSBM
  

[d, T] = size(X);

Q = eye(d);

end