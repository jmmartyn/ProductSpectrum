function y=Bvec2E_Min_Local(Bvec,M,rvec,n,Z,X,hz,hx,J,a,j)

%Converts Bvec to the contribution to energy that is affected by varying
%the unitary at site j (this unitary does not affect all contributions to
%the energy; it only affects those that are near site j, which is what we
%minimize over

%Constructs unitary from Bvec and updates the set of 2-qubit unitaries
V = makeUnitary(Bvec);
M{j,a} = V;
if j==n
  M{1,a} = V;
else
  M{j+1,a} = V;
end
  

low1 = j;       %lowest site where 1 site observables are affected by the unitary at site j
high1=j+1;      %highest site where 1 site observables are affected by the unitary at site j
if j==1
    low2=1;     %lowest site where 1 site observables are affected by the unitary at site j
    high2=j+1;  %highest site where 1 site observables are affected by the unitary at site j
elseif j==n-1
    low2=j-1;
    high2=n-1; 
else
    low2=j-1;
    high2=j+1; 
end


%Calculates the contribution to energy that is affected by the unitary at site j
E=0;
for i=low1:high1
    E=E-hz*LocalObs1(i, n, a, rvec, Z, Z, M)-hx*LocalObs1(i, n, a, rvec, X, Z, M);
end

for i=low2:high2
  E=E-J*LocalObs2(i, n, a, rvec, Z,Z, Z, M);
end 

y=E;
%y = (E+E')/2; %Removed to avoid strange divergences
