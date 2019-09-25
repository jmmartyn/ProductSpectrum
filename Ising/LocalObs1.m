function y=LocalObs1(j, n, current_layer, rvec, obs, Z, M)

%Calculates expectation value of obs_j (one qubit operator located at site j) using local techniques


X = [0 1; 1 0];


if 2*current_layer<n
  %Determines sites that are needed to compute expectation value locally
  if (j==1 || j==n) && mod(current_layer,2)==0
    current_layer=current_layer-1;
  end
  if mod(j,2)~= mod(current_layer,2) 
    q=j-1;
    left_site=q;
  else
    q=j; 
    left_site=q;
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
    s=j-1;
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
  %high=low+2*(a-1); %Highest left site of interest to local calculation
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
  %Ensures rho_bare_reduced is Hermitian and has unit trace
  rho_prod_reduced = (rho_prod_reduced+rho_prod_reduced')/2;
  rho_prod_reduced = rho_prod_reduced/trace(rho_prod_reduced);

  %Makes tensor product matrix of observable
  if j==low
    obs_matrix = obs;
    for i=low+1:high+1
        obs_matrix = kron(obs_matrix, eye(2,2));
    end
  else
    obs_matrix = eye(2,2);
    for i=low+1:j-1
      obs_matrix = kron(obs_matrix,eye(2,2));
    end
    obs_matrix =  kron(obs_matrix, obs);
    for i=j+1:high+1
      obs_matrix = kron(obs_matrix,eye(2,2));
    end
  end
  
  layer_matrices = cell(1,current_layer); %Product of unitary matrices at each layer
  
  layer_matrices{1} = cell2mat(M(low,1));
  for k=low+2:2:high
    layer_matrices{1} = kron(layer_matrices{1}, cell2mat(M(k,1)));
  end
  
  q=left_site;
  s=left_site;
  if current_layer==1
    layer_matrices{current_layer}=layer_matrices{1};
  elseif q==n-1 && current_layer==1
    layer_matrices{current_layer} = cell2mat(M(n-1,current_layer));  
  elseif q==1
    layer_matrices{current_layer} = cell2mat(M(1,current_layer));
    for i=3:high+1
        layer_matrices{current_layer} = kron(layer_matrices{current_layer}, eye(2,2));
    end
  else
    layer_matrices{current_layer} = eye(2,2);
    for i=low+1:left_site-1
      layer_matrices{current_layer} = kron(layer_matrices{current_layer},eye(2,2));
    end
    layer_matrices{current_layer} = kron(layer_matrices{current_layer}, cell2mat(M(left_site,current_layer)));
    for i=left_site+2:high+1
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
      for k=3:2:s
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


if 2*current_layer>=n 
  %Computes expectation value without local techniques if possible
  
  rho_prod = (eye(2,2) + rvec(1,1)*X+rvec(1,2)*Z)/2; %Product spectrum density matrix
  for i = 2:n
    rho_prod = kron(rho_prod, (eye(2,2)+rvec(i,1)*X+rvec(i,2)*Z)/2);
  end
  %Ensures rho_prod is Hermitian and has unit trace
  rho_prod = (rho_prod+rho_prod')/2;
  rho_prod = rho_prod/trace(rho_prod);
  
  %Constructs matrix obs_j
  if j==1
    obs_matrix = obs;
    for i=2:n
        obs_matrix = kron(obs_matrix,eye(2,2));
    end
  else
    obs_matrix = eye(2,2);
    for i=2:j-1
      obs_matrix = kron(obs_matrix,eye(2,2));
    end
    obs_matrix = kron(obs_matrix,obs);
    for i=j+1:n
      obs_matrix = kron(obs_matrix,eye(2,2));
    end
  end  
  
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
  
  %Product of layers of 2-qubit unitaries
  U=L{1};
  for i = 2:current_layer
    U = L{i}*U; %Be careful about commutativity here
  end    
    
  %Returns expectation value
  y = trace(obs_matrix*U*rho_prod*U');
  y = (y+y')/2;
end
