%Calculates product spectrum approximation to Ising thermal state

%Input parameters
tol = 5e-3;                                         %Tolerance of minimization
n=6;                                                %Number of sites; assumed to be even
l = 2;                                              %Number of layers
v = 2;                                              %Number of times each layer is varied over
num_vary = 2;                                       %Maximum number of rvec parameters varied over at once (2*num_vary real parameters are varied in total)
num_vary_times = 3;                                 %Number of times rvec is varied over
T = 1;                                              %Temperature
J=1;                                                %Coupling strength
hx=1.05;                                            %Magnetic field
hz=.5;                                              %Magnetic field
num_mins = (n/2*ceil(l/2)+(n/2-1)*floor(l/2))*v;    %Total number of minimizations


%Defines Pauli matrices
Z = [1 0; 0 -1];
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
I = [1 0; 0 1];
sigma = cell(1,4);
sigma{1} = I;
sigma{2} = X;
sigma{3} = Y;
sigma{4} = Z; 


%Finds exact density matrix and free energy for later comparison
%{
H = make_ising_hamiltonian(J,hx,hz,n,0);
[W,E]=eig(full(H));
rho_exact=expm(-E/T);
rho_exact=rho_exact/trace(rho_exact);
rho_exact = W*rho_exact*W';
F_exact = trace(H*rho_exact) + T*trace(rho_exact*logm(rho_exact));
%}


%Bloch vector parameters that specify the product state; we will record the set of parameters as each layer is varied over
rvecs = cell(1,l+1);


%Calculates the product state density matrix with minimum free energy
M = cell(n,1);          %2-qubit unitaries applied to product state, needed for variation over product state
for i=1:n
  M{i,1} = eye(4,4);    %Initially, the unitaries are the identity
end
current_layer=1;
rvec = zeros(n,2);      %n Bloch vector parameters used to determine density matrix
options = optimset('MaxFunEval', 100000000, 'TolFun', tol, 'TolX', tol, 'MaxIter',2*n*400, 'display','off');

%Minimizes free energy over product states
if num_vary >= n        %if (number of rvec parameters to be varied at once) > (total number of rvec parameters)
  %Minimizes free energy over all rvec parameters and renormalizes rvec if necessary
  [r0vec,F_min0]=fminunc(@(x)F_Local(T,M,x,n,Z,X,hz,hx,J,current_layer),rvec, options);
  for i=1:n
    if norm(r0vec(i,:))>1
      r0vec(i,:)=r0vec(i,:)/norm(r0vec(i,:));
    end
  end
  r1vec=r0vec;
  rvec=r1vec;
else %Otherwise, repeatedly minimize over only num_vary parameters at once; do this num_vary_times times 
  for current_variation=1:num_vary_times
    %rvec parameters from sites site1 to site2 are varied over 
    site1=1;
    site2=num_vary;
    for site_vary = 0:ceil(n/num_vary)-1
      r_vary_vec = zeros(site2-site1+1,2); %Parameters that are actually varied over
      for i=1:(site2-site1+1)
        r_vary_vec(i,1) = rvec(site_vary*num_vary+i,1);
        r_vary_vec(i,2) = rvec(site_vary*num_vary+i,2);
      end

      %Minimize over selected parameters and renormalizes rvec if necessary
      [r_vary0_vec,F_min0]=fminunc(@(x)F_Min_Local_Split(T,M,x,rvec, n,Z,X,hz,hx,J,current_layer,site1,site2),r_vary_vec, options);
      for i=1:length(r_vary0_vec)
        if norm(r_vary0_vec(i,:))>1
          r_vary0_vec(i,:)=r_vary0_vec(i,:)/norm(r_vary0_vec(i,:));
        end
      end
      r_vary1_vec=r_vary0_vec;
      for i=site1:site2
        rvec(i,1) = r_vary1_vec(i-site1+1,1);
        rvec(i,2) = r_vary1_vec(i-site1+1,2);
      end
      for i=1:n
        if norm(rvec(i,:))>1
          rvec(i,:)=rvec(i,:)/norm(rvec(i,:));
        end
      end

      %Increase values of site1 and site2 so that we vary over all rvec parameters
      site1=site1+num_vary; 
      if (site2+num_vary) <=n
        site2=site2+num_vary;
      else 
        site2=n;
      end
    end
  end
end

%Renormalizes rvec if necesssary
for i=1:n
  if norm(rvec(i,:))>1
   rvec(i,:)=rvec(i,:)/norm(rvec(i,:));
  end
end
rvecs{1} = rvec; %parameter values with no layers applied (l=0)



%Now apply unitaries and vary over both the unitaries and product state to minimize free energy
M = cell(n,l); %Array of unitary matrices applied to product spectrum
for i=1:n
    for layer=1:l
        M{i,layer} = eye(4,4);  %Initially, the unitaries are the identity 
    end
end

%Minimize energy over unitaries and product state
count=1; %Number of minimizations completed
for current_layer = 1:l
  start_pos = 1+((-1)^current_layer+1)/2; %Starting position; alternates between 1 and 2
  for b = 1:v
    for site = start_pos:2:n-1                         %site number
      Bvec=Unitary2Vec(M{site,current_layer});         %vector of parameters that parameterize the 2-qubit unitaries
      
      %Minimizes free energy over Bvec and determines the new 2-qubit unitary from Bvec
      [B0,E_min0]=fminunc(@(x)Bvec2E_Min_Local(x,M,rvec,n,Z,X,hz,hx,J,current_layer,site),Bvec, options); 
      Bvec = B0;
      V = makeUnitary(Bvec);
      M{site,current_layer} = V;
      if site==n
        M{1,current_layer} = V;
      else
        M{site+1,current_layer} = V;
      end
      
      %Minimizes free energy over rvec
      rvec = zeros(n,2); 
      if num_vary >= n    
        [r0vec,F_min0]=fminunc(@(x)F_Local(T,M,x,n,Z,X,hz,hx,J,current_layer),rvec, options);
        for i=1:n
          if norm(r0vec(i,:))>1
            r0vec(i,:)=r0vec(i,:)/norm(r0vec(i,:));
          end
        end
        rvec=r0vec;
      else
        for current_variation=1:num_vary_times
          site1=1;
          site2=num_vary;
          for site_vary = 0:ceil(n/num_vary)-1
            r_vary_vec = zeros(site2-site1+1,2);
            for i=1:(site2-site1+1)
              r_vary_vec(i,1) = rvec(site_vary*num_vary+i,1);
              r_vary_vec(i,2) = rvec(site_vary*num_vary+i,2);
            end
            [r_vary0_vec,F_min0]=fminunc(@(x)F_Min_Local_Split(T,M,x,rvec, n,Z,X,hz,hx,J,current_layer,site1,site2),r_vary_vec, options);
            for i=1:length(r_vary0_vec)
              if norm(r_vary0_vec(i,:))>1
                r_vary0_vec(i,:)=r_vary0_vec(i,:)/norm(r_vary0_vec(i,:));
              end
            end
            r_vary1_vec=r_vary0_vec;
            for i=site1:site2 
              rvec(i,1) = r_vary1_vec(i-site1+1,1);
              rvec(i,2) = r_vary1_vec(i-site1+1,2);
            end
            for i=1:n
              if norm(rvec(i,:))>1
                rvec(i,:)=rvec(i,:)/norm(rvec(i,:));
              end
            end
            site1=site1+num_vary; 
            if (site2+num_vary) <=n
              site2=site2+num_vary;
            else 
              site2=n;
            end
          end
        end
      end
      for i=1:n
        if norm(rvec(i,:))>1
          rvec(i,:)=rvec(i,:)/norm(rvec(i,:));
        end
      end
      disp(strcat(num2str(count),'/', num2str(num_mins),' minimizations completed.'));
      count=count+1;
    end
  end
  rvecs{current_layer+1} = rvec;
end

%Approximate minimum free energy
F_approx = F_Local(T,M,rvec,n,Z,X,hz,hx,J,current_layer);

%Analyzes approximate and exact thermal states
%{
total_matrix=make_total_matrix(M,l,n); 
rho_prod = (eye(2,2) + rvec(1,1)*X+rvec(1,2)*Z)/2;
for i = 2:n
  rho_prod = kron(rho_prod, (eye(2,2) + rvec(i,1)*X+rvec(i,2)*Z)/2);
end
%Ensures rho_bare is Hermitian and has unit trace
rho_prod = (rho_prod+rho_prod')/2;
rho_prod = rho_prod/trace(rho_prod);
rho_approx = total_matrix'*rho_prod*total_matrix;
rho_approx = rho_approx/trace(rho_approx); 
rho_approx = (rho_approx+rho_approx')/2;
F_dif = F_approx-F_exact;
F_frac_dif = abs(F_dif/F_exact);
%}


%Calculates approximate expectation values of spin
sigma_z_approx = zeros(n,1);
sigma_x_approx = zeros(n,1);
for  i=1:n
  sigma_z_approx(i) = LocalObs1(i, n, l, rvec, Z, Z, M);
  sigma_x_approx(i) = LocalObs1(i, n, l, rvec, X, Z, M);
end

%Plots spins
t = 1:n;

figure; scatter(t, sigma_z_approx, 'filled');
axis([1 n -.05 (max(sigma_z_approx)+.1)]);
leg=legend('Product Spectrum Z');
set(leg,'Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('$\sigma_z$','Interpreter','latex');
title('$\sigma_z$ comparison', 'Interpreter', 'latex')

figure; scatter(t, sigma_x_approx, 'filled');
axis([1 n -.05 (max(sigma_x_approx)+.1)]);
leg=legend('Product Spectrum X');
set(leg,'Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('$\sigma_x$','Interpreter','latex');
title('$\sigma_x$ comparison', 'Interpreter','latex')