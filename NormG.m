function g = NormG(V,R)
g = zeros(size(V,1),1);
for index = 1:size(V,1)
    v = V(index,:);
    g(index) = sum((v*R).^2);
end
end