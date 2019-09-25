function y=weyl(a,q)

%Constructs tensor product operator (X_i, Z_j, etc.) used in Ising Hamiltonian

omega=exp(2*pi*1i/q);

X=zeros(q,q);
for i=1:q
    for j=1:q
        if mod(i-j,q)==1
            X(i,j)=1;
        end
    end
end

Z=zeros(q,q);
for i=1:q
    Z(i,i)=omega^(i-1);
end

y=X^a(1)*Z^a(2);