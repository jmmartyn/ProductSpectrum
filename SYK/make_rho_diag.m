function y=make_rho_diag(r)
%Constructs product state from vector r of Bloch vector parameters

n=length(r);
rho_diag=sparse((1/2)*[1+r(1) 0; 0 1-r(1)]);
for i=2:n
    rho_diag=kron(sparse((1/2)*[1+r(i) 0; 0 1-r(i)]),rho_diag); %(I+r*Z)/2
end

y=rho_diag;

end