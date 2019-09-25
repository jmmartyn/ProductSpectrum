function y=Unitary2Vec(V)

%Converts unitary matrix to a vector of parameters (which parameterize the 
%unitary). This vector effectively parameterizes SU(4).

%Unitary is V = exp(1i*B); B is Hermitian
B = -1i*logm(V);
B = (B+B')/2;
Bvec = zeros(16,1);
count = 1;

%Diagonal terms
for i = 1:4
    Bvec(count)=B(i,i);
    count=count+1;
end

%Off-diagonal terms
for i = 1:4
  for j=i+1:4
    Bvec(count)=real(B(i,j));
    Bvec(count+6) = imag(B(i,j));
    count = count+1;
  end
end

y = Bvec;