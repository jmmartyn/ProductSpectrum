function y=F_Min_Local_Split(T,M,r_vary_vec,rvec, n,Z,X,hz,hx,J,current_layer,site1,site2)

%Calculates the contribution to the free energy from the sites site1 to
%site2

%Parameters that are varied (r_var_vec) and those that are not (rvec)
r_vec_parameters = zeros(n,2);
for i = 1:site1-1
    r_vec_parameters(i,1) = rvec(i,1);
    r_vec_parameters(i,2) = rvec(i,2);
end
for i=site1:site2
    r_vec_parameters(i,1) = r_vary_vec(i-site1+1,1);
    r_vec_parameters(i,2) = r_vary_vec(i-site1+1,2);
end
for i= site2+1:n
    r_vec_parameters(i,1) = rvec(i,1); 
    r_vec_parameters(i,2) = rvec(i,2); 
end

%Removes small imaginary components (rounding errors) and renormalizes rvec
%if necessary
r_vec_parameters = (r_vec_parameters+conj(r_vec_parameters))/2;
for i=1:n
  if norm(r_vec_parameters(i,:))>1
    r_vec_parameters(i,:)=r_vec_parameters(i,:)/norm(r_vec_parameters(i,:));
  end
end


%Determines the sites that are affected by the unitaries and product states from site1 through site2
if mod(site1,2)==0
    low=site1-1;
else
    low=site1;
end
if mod(site2,2)==0
    high=site2-1;
else
    high=site2;
end

for i=1:current_layer
    if low ==1
      break;
    else
      low=low-1;
    end
    if high ==n-1
      break;
    else
      high=high+1;
    end    
end

low1=low;
high1=high+1;

if low ==1
  low2=1;
else
  low2=low-1;
end  

if high == n-1
    high2 = n-1;
else
    high2=high+1;
end

%Entropy contribution
S = 0;
for i=site1:site2
  S = S -1/2*log((1+norm(r_vec_parameters(i,:)))^(1+norm(r_vec_parameters(i,:))))-1/2*log((1-norm(r_vec_parameters(i,:)))^(1-norm(r_vec_parameters(i,:))))+log(2);
end

%Energy contribution
E=0;
for i=low1:high1
    E=E-hz*LocalObs1(i, n, current_layer, r_vec_parameters, Z, Z, M)-hx*LocalObs1(i, n, current_layer, r_vec_parameters, X, Z, M);
end

for i=low2:high2
  E=E-J*LocalObs2(i, n, current_layer, r_vec_parameters, Z,Z, Z, M);
end 

%Returns free energy
y = E - T*S;
y=real(y)+imag(y);
%y = (y+y')/2; %Removed to avoid divergences