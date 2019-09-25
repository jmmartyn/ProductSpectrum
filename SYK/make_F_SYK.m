function F = make_F_SYK(r,E_list,T)

%Computes free energy of state with product state defined by r, and 2-qubit
%unitaries forming unitary
%Here, E_list = diag(unitary'*H*unitary)

%List of terms in product state
n=length(r);
rho_diag_list=(1/2)*[1+r(1),1-r(1)];
for ii=2:n
    rho_diag_list=kron((1/2)*[1+r(ii),1-r(ii)],rho_diag_list);
end
rho_diag_list=rho_diag_list';


%Entropy, unaffected by unitary
S=sum(-log(rho_diag_list.^rho_diag_list));

%Energy
E=sum(E_list.*rho_diag_list);

%Returns free energy
F=real(E-T*S)+abs(imag(E-T*S));


