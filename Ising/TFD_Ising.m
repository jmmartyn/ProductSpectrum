%Calculates TFD and expectation values (exact and approximate)
%This assumes r_vec and rho_exact are already known (which can be computed
%by executing the appropriate lines in product_spectrum_Ising)

obs1 = X;
j = 1; %Site that obs1 acts at; j<=k is asumed here
obs2 = X;
k = 2; %Site that obs2 acts at; j<=k is asumed here

%Constructs matrices corresponding to obs1_j and obs2_k
if j==1
    obs_matrix1 = obs1;
    for i=2:n
        obs_matrix1 = kron(obs_matrix1, eye(2,2));
    end
else
    obs_matrix1 = eye(2,2);
    for i=2:j-1
      obs_matrix1 = kron(obs_matrix1,eye(2,2));
    end
    obs_matrix1 =  kron(obs_matrix1, obs1);
    for i=j+1:n
      obs_matrix1 = kron(obs_matrix1,eye(2,2));
    end
end

if k==1
    obs_matrix2 = obs2;
    for i=2:n
        obs_matrix2 = kron(obs_matrix2, eye(2,2));
    end
else
    obs_matrix2 = eye(2,2);
    for i=2:k-1
      obs_matrix2 = kron(obs_matrix2, eye(2,2));
    end
    obs_matrix2 =  kron(obs_matrix2, obs2);
    for i=k+1:n
      obs_matrix2 = kron(obs_matrix2,eye(2,2));
    end
end

%calculates exact and approximate expectation values
TFD_exact = sqrtm(rho_exact); 
obs_exact = trace(TFD_exact'*obs_matrix1*TFD_exact*transpose(obs_matrix2));
obs_exact = real(obs_exact);
obs_approx = TFD_LocalObs2(j, k, n, l, rvec, obs1,obs2, Z, M);