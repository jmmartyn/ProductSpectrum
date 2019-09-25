function y=TFD_LocalObs2(j, k, n, a, rvec, obs1,obs2, Z, M)
%Computes <obs1_j*obs2_k> in TFD state using local techniques
%This assumes n is even
%This assumes j<=k

%
X = [0 1; 1 0];
if 2*a+2<n
  %Determines sites that are needed to compute expectation value locally
  if mod(j,2)~= mod(a,2) 
    if j==1
      q=2;
      left_site1=q;
      left_site2=q;
    elseif j==n-1
      q=j-1;
      left_site1=q;
      left_site2=q; 
    else
      q=j-1;
      left_site1=q;
      left_site2=j+1; 
    end
  else
    q=j; 
    left_site1=q;
    left_site2=q;
  end
  for i = 1:a-1
    if q==1
      break;
    else
      q=q-1;
    end
  end
  low=q; %Leftmost site of interest to local calculation
  
  if mod(k,2)~= mod(a,2) 
    if k==n-1
      s=n-1;
    elseif k==n
      s=n-1;  
    else
      s=k+1;
    end
  else
    s=k; 
  end
  for i = 1:a-1
    if s==n-1
      break;
    elseif s==n
      s=n-1; 
    else
      s=s+1;
    end
  end
  high=s; %Highest left site of interest to local calculation
  if left_site2==n+1
    left_site2=n-1;
  end

  num_sites = high-low+2; %Total number of sites of interest
  rvec_new = zeros(num_sites,2); %matrix of product spectrum parameters
  
  for i=1:num_sites
    rvec_new(i,1) = rvec(low+i-1,1);
    rvec_new(i,2) = rvec(low+i-1,2);
  end

  rho_bare_reduced = (eye(2,2) + rvec_new(1,1)*X+rvec_new(1,2)*Z)/2; %Product spectrum density matrix of reduced size
  for i = 2:num_sites
    rho_bare_reduced = kron(rho_bare_reduced, (eye(2,2) + rvec_new(i,1)*X+rvec_new(i,2)*Z)/2);
  end
  
  %Ensures rho_bare_reduced is Hermitian and has unit trace
  rho_bare_reduced = (rho_bare_reduced+rho_bare_reduced')/2;
  rho_bare_reduced = rho_bare_reduced/trace(rho_bare_reduced);

  %Makes tensor product matrix of observable
  if j==low
    obs_matrix1 = obs1;
    for i=low+1:high+1
        obs_matrix1 = kron(obs_matrix1, eye(2,2));
    end
  else
    obs_matrix1 = eye(2,2);
    for i=low+1:j-1
      obs_matrix1 = kron(obs_matrix1,eye(2,2));
    end
    obs_matrix1 = kron(obs_matrix1, obs1);
    for i=j+1:high+1
      obs_matrix1 = kron(obs_matrix1,eye(2,2));
    end
  end

  if k==low
    obs_matrix2 = obs2;
    for i=low+1:high+1
        obs_matrix2 = kron(obs_matrix2, eye(2,2));
    end
  else
    obs_matrix2 = eye(2,2);
    for i=low+1:k-1
      obs_matrix2 = kron(obs_matrix2,eye(2,2));
    end
    obs_matrix2 = kron(obs_matrix2, obs2);
    for i=k+1:high+1
      obs_matrix2 = kron(obs_matrix2,eye(2,2));
    end
  end

  %Computes full unitary matrices at each layer
  layer_matrices = cell(1,a);
  layer_matrices{1} = cell2mat(M(low,1));
  for k=low+2:2:high
    layer_matrices{1} = kron(layer_matrices{1}, cell2mat(M(k,1)));
  end
  
  q=left_site1;
  s=left_site2;

  if a==1
    layer_matrices{a}=layer_matrices{1};
  elseif q==n-1 && a==1
    layer_matrices{a} = cell2mat(M(n-1,a));
  elseif q==n && mod(a,2)==0
    sz = size(layer_matrices{1});
    layer_matrices{a} = eye(sz(1),sz(2)); 
  elseif q==1
    layer_matrices{a} = cell2mat(M(1,a));
    for k=q+2:2:s
        layer_matrices{a} = kron(layer_matrices{a}, cell2mat(M(k,a)));
    end
    for k=s+2:high+1
        layer_matrices{a} = kron(layer_matrices{a}, eye(2,2));
    end
  else
    layer_matrices{a} = eye(2,2);
    for k=low+1:q-1
      layer_matrices{a} = kron(layer_matrices{a},eye(2,2));
    end
    for k=q:2:s
      layer_matrices{a} = kron(layer_matrices{a}, cell2mat(M(k,a)));
    end    
    for k=s+2:high+1
      layer_matrices{a} = kron(layer_matrices{a},eye(2,2));
    end
  end
  if q==1 && mod(a,2)==1
      q=2;
  elseif q==1 && mod(a,2)==0
      q=1;
  else
      q=q-1;
  end
  if s==n-1 && mod(a,2)==1
      s=n-2;
  elseif s==n-1 && mod(a,2)==0
      s=n-1;
  else
      s=s+1;
  end
  
  for i=a-1:-1:2
    if q==1
      layer_matrices{i} = cell2mat(M(1,i));
      for k=q+2:2:s
        layer_matrices{i}=kron(layer_matrices{i}, cell2mat(M(k,i)));    
      end
      for k=s+2:high+1
        layer_matrices{i} = kron(layer_matrices{i}, eye(2,2));
      end
    else
      layer_matrices{i} = eye(2,2);
      for k=low+1:q-1
        layer_matrices{i} = kron(layer_matrices{i},eye(2,2));
      end
      for k=q:2:s
        layer_matrices{i} =  kron(layer_matrices{i}, cell2mat(M(k,i)));
      end
      for k=s+2:high+1
        layer_matrices{i} = kron(layer_matrices{i},eye(2,2));
      end
    end
    if q==1 && mod(i,2)==1
      q=2;
    elseif q==1 && mod(i,2)==0
      q=1;
    else
      q=q-1;
    end
    if s==n-1 && mod(i,2)==1
      s=n-2;
    elseif s==n-1 && mod(i,2)==0
      s=n-1;
    else
      s=s+1;
    end   
  end
  
  U =layer_matrices{1};
  for i=2:a
    U = layer_matrices{i}*U; %Be careful about commutativity here
  end
  
  TFD_bare_reduced = sqrtm(rho_bare_reduced);
  TFD_reduced = U*TFD_bare_reduced*U';
  
  %Returns expectation value
  y = trace(TFD_reduced*obs_matrix1*TFD_reduced*transpose(obs_matrix2));
  y=(y+y')/2;
end




if (2*a+2) >= n 
  %Computes expectation value without local techniques if possible
  
  rho_bare = (eye(2,2) + rvec(1,1)*X+rvec(1,2)*Z)/2; %Product spectrum density matrix
  for i = 2:n
    rho_bare = kron(rho_bare, (eye(2,2)+rvec(i,1)*X+rvec(i,2)*Z)/2);
  end
  
  %Ensures rho_bare is Hermitian and has unit trace
  rho_bare = (rho_bare+rho_bare')/2;
  rho_bare = rho_bare/trace(rho_bare);
  
  %Constructs matrix obs1_j*obs2_(j+1)
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
      obs_matrix2 = kron(obs_matrix2,eye(2,2));
    end
    obs_matrix2 =  kron(obs_matrix2, obs2);
    for i=k+1:n
      obs_matrix2 = kron(obs_matrix2,eye(2,2));
    end
  end

  %obs_matrix = obs_matrix1*obs_matrix2;
  
  %Layer of 2 qubit unitaries
  L = cell(1,a);
  for i=1:a
    p = 1+((-1)^i+1)/2;
    L{i} = M{p,i};
    for c = (p+2):2:n-1
        L{i} = kron(L{i}, M{c,i});
    end
    if mod(p,2)==0
        L{i} = kron(eye(2,2), kron(L{i}, eye(2,2)));
    end
  end
  
  %Constructs product of layer_matrices
  U=L{1};
  for i = 2:a
    U = L{i}*U; %Be careful about commutativity here
  end    
  
  %Constructs TFD
  TFD_bare = sqrtm(rho_bare);
  TFD = U*TFD_bare*U';
  
  %Returns expectation value
  y = trace(TFD*obs_matrix1*TFD*transpose(obs_matrix2));
  y = (y+y')/2;
end
