function y=globalweyl(r,a,n,q)

%Constructs tensor product operator (X_i, Z_j, etc.) used in Ising Hamiltonian
if r==1
    y=weyl(a,q);
else
    y=eye(q);
end

for i=2:n
    if i==r
        y=kron(y,weyl(a,q));
    else
        y=kron(y,eye(q));
    end
end

y=sparse(y);

end

    

