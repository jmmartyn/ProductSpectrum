function y=makeUnitary(Bvec)

%Converts Bvec to 2-qubit unitary matrix

B = zeros(4,4); %B will be Hermitian
count = 1;

%Diagonal terms
for i = 1:4
    B(i,i) = Bvec(count);
    count = count+1; 
end

%Off-diagonal terms
for i = 1:4
  for j=i+1:4
    B(i,j) = Bvec(count)+1i*Bvec(count+6);
    count = count+1;
  end
end
for i = 1:4
  for j=i:4
    B(j,i) = B(i,j)';
  end
end

%Ensures B is Hermitian
B = (B+B')/2;

%Returns the unitary
y = expm(1i*B);



