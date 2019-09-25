function y=F_Local(T,M,rvec,n,Z,X,hz,hx,J,a)

%Converts parameters, rvec, and unitaries (in M) to free energy

%Removes small imaginary components (rounding errors) and renormalizes rvec
%if necessary
rvec = (rvec+conj(rvec))/2; 
for i=1:n
  if norm(rvec(i,:))>1
    rvec(i,:)=rvec(i,:)/norm(rvec(i,:));
  end
end


%Entropy of product state
S = 0;
for i=1:n
  S = S -1/2*log((1+norm(rvec(i,:)))^(1+norm(rvec(i,:))))-1/2*log((1-norm(rvec(i,:)))^(1-norm(rvec(i,:))))+log(2);
end

%Energy of state
E=0;
for j=1:n-1
    E=E-J*LocalObs2(j, n, a, rvec, Z,Z, Z, M);
end
for j=1:n
    E=E-hz*LocalObs1(j, n, a, rvec, Z, Z, M)-hx*LocalObs1(j, n, a, rvec, X, Z, M);
end

%Returns free energy
y = E - T*S;
%y = (y+y')/2; %Removed to avoid divergences