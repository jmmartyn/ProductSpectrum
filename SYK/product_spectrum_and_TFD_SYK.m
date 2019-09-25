%Calculates product spectrum approximation to SYK thermal state and SYK TFD

%Input parameters
n=4;                    %Number of spins, assumed to be even
D=2^n;                  %Dimension
m=2*n;                  %Number of spins in TFD
T=1;                    %Temperature
J=1;                    %Coupling
num_layers = 6;         %Number of layers
rep_tot = 10;           %Number of repetitions (sweeps over circuit)
i_tot = 4;              %Number of iterations over a single unitary
rand_multiple = .3;     %Random number used to determine starting unitary for minimizations
options = optimoptions('fmincon','MaxIterations',1e10,'TolFun', 1e-8, ...
    'TolX', 1e-8,'Display','off','MaxFunctionEvaluations',1e10); %Minimization options
sw0=1;                  %Sets up new Hamiltonian if 1
sw1=1;                  %Sets up new permutations if 1


%Sets up Hamiltonian
if sw0==1 % if sw0=0 then the Hamiltonian etc in memory is used
  %Sets up parallel pool for parellel computations
  p = gcp('nocreate'); %If no pool, do not create new one.
  if isempty(p)
    parpool;
  end
  
  %Computes SYK Hamiltonian and additional variables
  [X, Y, Z, chi, H]=make_SYK_Hamiltonian_and_vars(n,m,D);  
  
  %Exact solution by diagonalization
  [B, E]=eig(full(H));
  E_list=zeros(D,1);
  for k=1:D
    E_list(k)=E(k,k);
  end
end

%Constructs exact thermal state
low = min(E_list);
E_list = E_list-low;
rho_list=exp(-E_list/T);
rho_list=rho_list/sum(rho_list);
rho_diag = diag(rho_list);
H_new = H-E(2^n,2^n)*eye(2^n,2^n);  %Hamiltonian minus its largest eigenvalue


%Parameters in product spectrum state
r_list = ones(n,1);         %Bloch vector parameters
M = cell(n/2,num_layers);   %Set of 2-qubit unitaries
L = cell(1,num_layers);     %Layers of 2-qubit unitaries
for j=1:num_layers
  L{j} = eye(2^n,2^n); 
  for i=1:n/2
    M{i,j} = eye(4,4); 
  end
end

%Performs minimization over free energy
perm=PermutationOperator([2,2^(n-2),2],[3,1,2]);    %Operator used to shift operators and state by 2 spins   
for rep = 1:rep_tot
  rho_diag_approx = make_rho_diag(r_list);          %Makes product state
  for layer = 1:num_layers
    if mod(layer,2) ==1                             %if layer is odd
      layers_matrix_outer = eye(2^n,2^n);           %2-qubit unitaries above the current layer
      layers_matrix_inner = eye(2^n,2^n);           %2-qubit unitaries below the current layer
      for a=layer+1:num_layers
        layers_matrix_outer=L{a}*layers_matrix_outer;
      end
      for a=1:layer-1
        layers_matrix_inner=L{a}*layers_matrix_inner;
      end 
      
      %iteratively minimizes over unitaries using the iterative algorithm from arxiv:0707.1454
      for position = 1:2:n                  %minimize over 2-qubit unitaries at each position
        if rep ==0
          U_tmp = rand(4,4)+1i*rand(4,4);   %temporary 2-qubit unitary used for minimization
        else
          U_tmp = eye(4,4)+rand_multiple*(rand(4,4)+1i*rand(4,4));
        end
          
        %Matrices used to determine energy
        mat1 = L{layer}*(layers_matrix_inner*(rho_diag_approx*(layers_matrix_inner'*L{layer}')));
        mat2=layers_matrix_outer'*(H_new*layers_matrix_outer);
          
        %Performs iterative algorithm i_tot times
        for i=1:i_tot
          U_tmp_full = sparse(kron(eye(2^(position-1),2^(position-1)),kron(U_tmp,eye(2^(n-position-1),2^(n-position-1)))));
          Y_U = mat1*(U_tmp_full'*mat2);          %Environment of 2-qubit unitary under focus
          Y_U_red = PartialTrace(Y_U,[1,3],[2^(position-1),4,2^(n-position-1)]);
          [V,S,W]=svd(Y_U_red);
          U_tmp = -1*(W*V');                      %Unitary that minimizes energy
          U_tmp = U_tmp*(det(U_tmp))^(-2^(-n));   %Ensures det(U)=1
        end
          
        %Updates unitaries and layer of unitaries
        M{(position-1)/2+1,layer} = U_tmp*M{(position-1)/2+1,layer};
        U_tmp_full = sparse(kron(eye(2^(position-1),2^(position-1)),kron(U_tmp,eye(2^(n-position-1),2^(n-position-1)))));
        L{layer} = U_tmp_full*L{layer};
      end
      
      %Updates unitaries applied to product state
      L{layer} = M{1,layer};
      for a=2:n/2
        L{layer} = kron(L{layer},M{a,layer});
      end
      unitary = eye(2^n,2^n);
      for a=1:num_layers
        unitary=L{a}*unitary;
      end
      
      %rho_approx = unitary*(rho_diag_approx*unitary');
      %trace(rho_approx*H_new)
      
      %Minimizes free energy over Bloch vector parameters r_list
      Htmp=unitary'*H*unitary;
      Etmp=zeros(2^n,1);
      for a=1:2^n
        Etmp(a)=Htmp(a,a);
      end
      r_list_tmp=r_list;
      Ftmp=@(r)make_F_SYK(r, Etmp,T); 
      [r_list, F]=fmincon(Ftmp,r_list_tmp,zeros(n,n),zeros(n,1),zeros(n,n),zeros(n,1),-1*ones(n,1),ones(n,1),[],options);
      rho_diag_approx = make_rho_diag(r_list);
      rho_approx = unitary*(rho_diag_approx*unitary');
      %trace(rho_approx*H_new)
      %F
      
    else    %if layer is even, shift the unitaries by one spin, then apply the same minimization procedure
      layers_matrix_outer = eye(2^n,2^n);           %2-qubit unitaries above the current layer
      layers_matrix_inner = eye(2^n,2^n);           %2-qubit unitaries below the current layer
      for a=layer+1:num_layers
        layers_matrix_outer=L{a}*layers_matrix_outer;
      end
      for a=1:layer-1
        layers_matrix_inner=L{a}*layers_matrix_inner;
      end 
      
      %iteratively minimizes over unitaries using the iterative algorithm from arxiv:0707.1454
      for position = 1:2:n                  %minimize over 2-qubit unitaries at each position
        if rep ==0
          U_tmp = rand(4,4)+1i*rand(4,4);   %temporary 2-qubit unitary used for minimization
        else
          U_tmp = eye(4,4)+rand_multiple*(rand(4,4)+1i*rand(4,4));
        end
        
        %Matrices used to determine energy
        mat1=perm*(L{layer}*(layers_matrix_inner*(rho_diag_approx*(layers_matrix_inner'*(L{layer}'*perm')))));
        mat2=perm*(layers_matrix_outer'*(H_new*(layers_matrix_outer*perm')));
        
        %Performs iterative algorithm i_tot times
        for i=1:i_tot
          U_tmp_full = sparse(kron(eye(2^(position-1),2^(position-1)),kron(U_tmp,eye(2^(n-position-1),2^(n-position-1)))));
          Y_U = mat1*(U_tmp_full'*mat2);            %Environment of 2-qubit unitary under focus
          Y_U_red = PartialTrace(Y_U,[1,3],[2^(position-1),4,2^(n-position-1)]);
          [V,S,W]=svd(Y_U_red);
          U_tmp = -1*(W*V');                        %Unitary that minimizes energy
          U_tmp = U_tmp*(det(U_tmp))^(-2^(-n));     %Ensures det(U)=1
        end
        
        %Updates unitaries and layer of unitaries
        M{(position-1)/2+1,layer} = U_tmp*M{(position-1)/2+1,layer};
        U_tmp_full = sparse(kron(eye(2^(position-1),2^(position-1)),kron(U_tmp,eye(2^(n-position-1),2^(n-position-1)))));
        L{layer} = perm'*(U_tmp_full*(perm*L{layer}));
      end
      
      %Updates unitaries applied to product state
      L{layer} = M{1,layer};
      for a=2:n/2
        L{layer} = kron(L{layer},M{a,layer});
      end
      L{layer} = perm'*(L{layer}*perm); 
      unitary = eye(2^n,2^n);
      for a=1:num_layers
        unitary=L{a}*unitary;
      end
      
      %rho_approx = unitary*(rho_diag_approx*unitary');
      %trace(rho_approx*H_new)
      
      %Minimizes free energy over Bloch vector parameters r_list
      Htmp=unitary'*H*unitary;
      Etmp=zeros(2^n,1);
      for a=1:2^n
        Etmp(a)=Htmp(a,a);
      end
      r_list_tmp=r_list;
      Ftmp=@(r)make_F_SYK(r, Etmp,T); 
      [r_list, F]=fmincon(Ftmp,r_list_tmp,zeros(n,n),zeros(n,1),zeros(n,n),zeros(n,1),-1*ones(n,1),ones(n,1),[],options);
      rho_diag_approx = make_rho_diag(r_list);
      rho_approx = unitary*(rho_diag_approx*unitary');
      
      %trace(rho_approx*H_new)
      %F
    end
  end
  disp(strcat('repetition: ', num2str(rep), '/', num2str(rep_tot)));
end 

%Exact and approximate free energy and density matrices
rho_exact = B*rho_diag*B';
rhosqrt_exact=B*diag(sqrt(rho_list))*B';
rhosqrt_approx=sqrtm(rho_approx);
F_approx=F; 
F_exact=sum(E_list.*rho_list)-T*sum(-log(rho_list.^rho_list));
FFracDif = abs((F_exact-F_approx)/F_exact);


%Time setup
t_min=0;
t_max=10;
t_num=4*(t_max-t_min);
t_list=linspace(t_min,t_max,t_num);

%One sided correlators and tfd correlators
chi2_exact=zeros(2*n,t_num);
chi2_approx=zeros(2*n,t_num);
chi2tfd_exact=zeros(2*n,t_num);
chi2tfd_approx=zeros(2*n,t_num);
for nt=1:t_num
    utmp=expm(-1i*H*t_list(nt));
    parfor i=1:2*n
        chi2_exact(i,nt)=trace(rho_exact*(utmp'*(chi{i}*(utmp*chi{i}))));
        chi2_approx(i,nt)=trace(rho_approx*(utmp'*(chi{i}*(utmp*chi{i}))));
        chi2tfd_exact(i,nt)=trace(rhosqrt_exact*(utmp'*(chi{i}*(utmp*(rhosqrt_exact*chi{i})))));
        chi2tfd_approx(i,nt)=trace(rhosqrt_approx*(utmp'*(chi{i}*(utmp*(rhosqrt_approx*chi{i})))));
    end
end

%Entropy
S_exact=zeros(n,t_num);
S_approx=zeros(n,t_num);
for nt=1:t_num
    utmp=expm(-1i*H*t_list(nt));
    for i=1:n
        if nt==1
            S_exact(i,nt)=VNent(TrX(rho_exact,[2],[2^i 2^(n-i)]));
        else
            S_exact(i,nt)=S_exact(i,1);
        end
        rhotmp=utmp*rho_approx*utmp';
        S_approx(i,nt)=VNent(TrX(rhotmp,[2],[2^i 2^(n-i)]));
    end
end

%Makes plots
col1=hsv(2*n);
figure();
hold on
for i=1:2*n
    plot(t_list,real(chi2tfd_exact(i,:)),'color',col1(i,:),'LineWidth',2);
    scatter(t_list,real(chi2tfd_approx(i,:)),[],col1(i,:),'LineWidth',1.3);
end
hold off
ylabel('$\langle \chi_\alpha(t) \otimes  \chi_\alpha^* \rangle$','Interpreter','latex','FontSize',24);
xlabel('$t$','Interpreter','latex','FontSize',24);
legend

figure();
hold on
for i=1:2*n
    plot(t_list,real(chi2_exact(i,:)),'color',col1(i,:),'LineWidth',2);
    scatter(t_list,real(chi2_approx(i,:)),[],col1(i,:),'LineWidth',1.3);
end
hold off
ylabel('$\langle \chi_\alpha(t) \chi_\alpha \rangle$','Interpreter','latex','FontSize',24);
xlabel('$t$','Interpreter','latex','FontSize',24);
legend

col2=hsv(n);
figure();
hold on
for i=n:-1:1
    plot(t_list,real(S_exact(i,:)),'color',col2(i,:),'LineWidth',2);
    scatter(t_list,real(S_approx(i,:)),[],col2(i,:),'LineWidth',1.3);
end
hold off
ylabel('$S$','Interpreter','latex','FontSize',24);
xlabel('$t$','Interpreter','latex','FontSize',24);
legend