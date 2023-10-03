function [ L ] = mblockchol( M, n, m )
%MBLOCKCHOL Block Cholesky M = L*L'
% n = size of each block
% m = number of blocks
% Alexander Mamonov, 2015
%==========================================================================

L = zeros(n*m, n*m);
ind = @(j) (1:n) + n*(j-1);

for k = 1:m
    msum = zeros(n, n);
    
	for j = 1:(k-1)
        msum = msum + L(ind(k), ind(j)) * L(ind(k), ind(j))';
    end
    
	L(ind(k), ind(k)) = sqrtm( M(ind(k), ind(k)) - msum );
    
    %L(ind(k), ind(k)) = real(sqrtm( M(ind(k), ind(k)) - msum ));
        
    for i = (k+1):m
        msum = zeros(n, n);
        
        for j = 1:(k-1)
            msum = msum + L(ind(i), ind(j)) * L(ind(k), ind(j))';
        end
        
        L(ind(i), ind(k)) =  ...
            (M(ind(i), ind(k)) - msum) / L(ind(k), ind(k));
        
        %L(ind(i), ind(k)) =  ...
        %    (M(ind(i), ind(k)) - msum) * pinv(L(ind(k), ind(k)));
        
    end
end

end

