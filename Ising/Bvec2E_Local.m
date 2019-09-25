function y=Bvec2E_Local(Bvec,M,rvec,n,Z,X,hz,hx,J,a,j)

%Converts Bvec to expected energy E

%Constructs unitary from Bvec and updates the set of 2-qubit unitaries
V = makeUnitary(Bvec);
M{j,a} = V;
if j==n
  M{1,a} = V;
else
  M{j+1,a} = V;
end

%Calculates energy
E=0;
for j=1:n-1
    E=E-J*LocalObs2(j, n, a, rvec, Z,Z, Z, M);
end
for j=1:n
    E=E-hz*LocalObs1(j, n, a, rvec, Z, Z, M)-hx*LocalObs1(j, n, a, rvec, X, Z, M);
end

y=E;
%y = (E+E')/2; %Removed to avoid strange divergences
