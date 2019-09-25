function y=s_weyl(a)

%Constructs sparse tensor product operator (X_i, Z_j, etc.) used in SYK Hamiltonian

X=[0 1; 1 0];
Z=[1 0; 0 -1];

y=X^a(1)*Z^a(2);

y=sparse(y);

end