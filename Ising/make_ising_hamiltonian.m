function y=make_ising_hamiltonian(J,hx,hz,n,BC)

%Make certain special Pauli operators; assumes n is even
Z=cell(1,n);
for k=1:n
    Z{k}=globalweyl(k,[0 1],n,2);
end
X=cell(1,n);
for k=1:n
    X{k}=globalweyl(k,[1 0],n,2);
end

%Make the Ising Hamiltonian
y=zeros(2^n,2^n);

for i=1:n-1
    y=y-J*Z{i}*Z{i+1};
end

if BC ==1 %Periodic boundary conditions
  y = y-J*Z{n}*Z{1};
end

for i=1:n
    y=y-hx*globalweyl(i,[1 0],n,2);
    y=y-hz*Z{i};
end

% remove roundoff errors (make hermitian)
y=(y+y')/2;

end