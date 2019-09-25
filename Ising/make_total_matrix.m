function y=make_total_matrix(M,l,n)

%Converts set of 2-qubit unitaries in M to their product over all sites and
%layers (this is the total unitary U that is applied to the product state)


%Product of matrices at each layer
L = cell(1,l);

for i=1:l
    p = 1+((-1)^i+1)/2;
    layerL{i} =cell2mat(M(p,i));
    for j = (p+2):2:n-1
        layerL{i} = kron(layerL{i}, cell2mat(M(j,i)));
    end
    if mod(p,2)==0
        layerL{i} = kron(kron(eye(2,2),layerL{i}),eye(2,2));
    end
end

y=layerL{1};
for i = 2:l
    y = layerL{i}*y;
end
