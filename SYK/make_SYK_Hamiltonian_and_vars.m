function [X, Y, Z, chi, H]=make_SYK_Hamiltonian_and_vars(n,m,D)
%Computes SYK Hamiltonian (H), Majorana fermions (chi), and Paulis (X, Y, Z)

% Pauli matrices
X=cell(n,1);
Y=cell(n,1);
Z=cell(n,1);
for r=1:n
    X{r}=s_global_weyl(r,[1 0],n);
    Y{r}=1i*s_global_weyl(r,[1 1],n);
    Z{r}=s_global_weyl(r,[0 1],n);
end

%Array of Majorana fermions, using Jordan-Wigner transformation
chi=cell(2*n,1);
for r=1:n
    chi{2*(r-1)+1}=speye(D);
    chi{2*(r-1)+2}=speye(D);
    for rp=1:r-1
        chi{2*(r-1)+1}=Z{rp}*chi{2*(r-1)+1};
        chi{2*(r-1)+2}=Z{rp}*chi{2*(r-1)+2};
    end
    for a=1:2
        if a==1
            chi{2*(r-1)+a}=X{r}*chi{2*(r-1)+a};
        else
            chi{2*(r-1)+a}=Y{r}*chi{2*(r-1)+a};
        end
    end
end
%Normalization of Majoranas
for r=1:2*n
    chi{r}=chi{r}/sqrt(2);
end

%Constructs terms in Hamiltonian
term = cell(2*n,1);
parfor i = 1:2*n
    term{i,1}=sparse(zeros(2^n,2^n));
end
J=1;
H=sparse(D,D);
parfor i1=1:2*n
    for i2=1:i1-1
        for i3=1:i2-1
            for i4=1:i3-1
                term{i1,1} = term{i1,1}+sparse(normrnd(0,J*sqrt(6/m^3))*chi{i4}*chi{i3}*chi{i2}*chi{i1});
            end
        end
    end
end
%Construct Hamiltonian
parfor i1=1:2*n
    H = H+term{i1,1};
end
H=sparse(H);



end