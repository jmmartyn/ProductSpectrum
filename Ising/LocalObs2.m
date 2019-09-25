function y=LocalObs2(j, n, current_layer, rvec, obs1,obs2, Z, M)
%Computes <obs1_j*obs2_(j+1)> using local techniques

X = [0 1; 1 0];

if 2*current_layer+2<n
  %Determines sites that are needed to compute expectation value locally
  if mod(j,2)~= mod(current_layer,2) 
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
  for i = 1:current_layer-1
    if q==1
      break;
    else
      q=q-1;
    end
  end
  low=q; %Leftmost site of interest to local calculation
  
  if mod(j,2)~= mod(current_layer,2) 
    if j==n-1
      s=n-1;
    else
      s=j+1;
    end
  else
    s=j; 
  end
  for i = 1:current_layer-1
    if s==n-1
      break;
    else
      s=s+1;
    end
  end
  high=s; %Highest left site of interest to local calculation
  
  num_sites = high-low+2; %Total number of sites of interest
  rvec_new = zeros(num_sites,2); %matrix of product spectrum parameters
  
  for i=1:num_sites
    rvec_new(i,1) = rvec(low+i-1,1);
    rvec_new(i,2) = rvec(low+i-1,2);
  end

  rho_prod_reduced = (eye(2,2) + rvec_new(1,1)*X+rvec_new(1,2)*Z)/2; %Product spectrum density matrix of reduced size
  for i = 2:num_sites
    rho_prod_reduced = kron(rho_prod_reduced, (eye(2,2) + rvec_new(i,1)*X+rvec_new(i,2)*Z)/2);
  end
  
  %Ensures rho_prod_reduced is Hermitian and has unit trace
  rho_prod_reduced = (rho_prod_reduced+rho_prod_reduced')/2;
  rho_prod_reduced = rho_prod_reduced/trace(rho_prod_reduced);

  %Makes tensor product matrix of observables
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

  obs_matrix2 = eye(2,2);
  for i=low+1:j
    obs_matrix2 = kron(obs_matrix2,eye(2,2));
  end
  obs_matrix2 = kron(obs_matrix2, obs2);
  for i=j+2:high+1
    obs_matrix2 = kron(obs_matrix2,eye(2,2));
  end
  obs_matrix = obs_matrix1*obs_matrix2;

  %Computes full unitary matrices at each layer
  layer_matrices = cell(1,current_layer); 
  layer_matrices{1} = cell2mat(M(low,1));
  for k=low+2:2:high
    layer_matrices{1} = kron(layer_matrices{1}, cell2mat(M(k,1)));
  end
  
  q=left_site1;
  s=left_site2;

  if current_layer==1
    layer_matrices{current_layer}=layer_matrices{1};
  elseif q==n-1 && current_layer==1
    layer_matrices{current_layer} = cell2mat(M(n-1,current_layer));
  elseif q==1
    layer_matrices{current_layer} = cell2mat(M(1,current_layer));
    for k=q+2:2:s
        layer_matrices{current_layer} = kron(layer_matrices{current_layer}, cell2mat(M(k,current_layer)));
    end
    for k=s+2:high+1
        layer_matrices{current_layer} = kron(layer_matrices{current_layer}, eye(2,2));
    end
  else
    layer_matrices{current_layer} = eye(2,2);
    for k=low+1:q-1
      layer_matrices{current_layer} = kron(layer_matrices{current_layer},eye(2,2));
    end
    for k=q:2:s
      layer_matrices{current_layer} = kron(layer_matrices{current_layer}, cell2mat(M(k,current_layer)));
    end    
    for k=s+2:high+1
      layer_matrices{current_layer} = kron(layer_matrices{current_layer},eye(2,2));
    end
  end
  if q==1 && mod(current_layer,2)==1
      q=2;
  elseif q==1 && mod(current_layer,2)==0
      q=1;
  else
      q=q-1;
  end
  if s==n-1 && mod(current_layer,2)==1
      s=n-2;
  elseif s==n-1 && mod(current_layer,2)==0
      s=n-1;
  else
      s=s+1;
  end
  
  for i=current_layer-1:-1:2
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

  %Constructs product of layer_matrices
  U =layer_matrices{1};
  for i=2:current_layer
    U = layer_matrices{i}*U; %Be careful about commutativity here
  end
  
  %Returns expectation value
  y = trace(obs_matrix*U*rho_prod_reduced*U');
  y=(y+y')/2;
end




if (2*current_layer+2) >= n
  %Computes expectation value without local techniques if possible
  
  rho_prod = (eye(2,2) + rvec(1,1)*X+rvec(1,2)*Z)/2; %Product spectrum density matrix
  for i = 2:n
    rho_prod = kron(rho_prod, (eye(2,2)+rvec(i,1)*X+rvec(i,2)*Z)/2);
  end
  %Ensures rho_prod is Hermitian and has unit trace
  rho_prod = (rho_prod+rho_prod')/2;
  rho_prod = rho_prod/trace(rho_prod);
  
  %Constructs matrix obs1_j*obs2_(j+1)
  low=1;
  high=n-1;
  if j==1
    obs_matrix1 = obs1;
    for i=2:high+1
        obs_matrix1 = kron(obs_matrix1, eye(2,2));
    end
  else
    obs_matrix1 = eye(2,2);
    for i=low+1:j-1
      obs_matrix1 = kron(obs_matrix1,eye(2,2));
    end
    obs_matrix1 =  kron(obs_matrix1, obs1);
    for i=j+1:high+1
      obs_matrix1 = kron(obs_matrix1,eye(2,2));
    end
  end

  obs_matrix2 = eye(2,2);
  for i=low+1:j
    obs_matrix2 = kron(obs_matrix2,eye(2,2));
  end
  obs_matrix2 =  kron(obs_matrix2, obs2);
  for i=j+2:high+1
    obs_matrix2 = kron(obs_matrix2,eye(2,2));
  end
  obs_matrix = obs_matrix1*obs_matrix2;
  
  %Layer of 2 qubit unitaries
  L = cell(1,current_layer);
  for i=1:current_layer
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
  for i = 2:current_layer
    U = L{i}*U; %Be careful about commutativity here
  end    
  
  %Returns expectation value
  y = trace(obs_matrix*U*rho_prod*U');
  y = (y+y')/2;
end
