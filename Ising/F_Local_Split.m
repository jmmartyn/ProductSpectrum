function y=F_Local_Split(T,M,r_vary_vec,rvec, n,Z,X,hz,hx,J,a,site1,site2)

%Calculates free energy from rvec, while varying only parameters from rvec(site1) to rvec(site2)

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

%Returns free energy
y= F_Local(T,M,r_vec_parameters,n,Z,X,hz,hx,J,a);
y=real(y)+imag(y); %ensures stability of algorithm