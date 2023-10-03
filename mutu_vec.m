function UTU = mutu_vec(gcos, n, nt)
%MUTPU assemple U'*U matrix from the data
% Alexander Mamonov, 2015
%==========================================================================
% Adapted to avoind cell structures and enable GPU

m1 = nt / 2;

UTU = zeros(n*m1, n*m1);

ind = @(j) n*(j-1) + (1:n);

for i = 1:m1
    for j = 1:m1
        
        UTU(ind(i), ind(j)) = ...
            (gcos(:,:,j + i - 1) + gcos(:,:, abs( j - i ) + 1));
%           reshape((gcos(:,j + i - 1) + gcos(:, abs( j - i ) + 1)),[n n]);
            
      
    end
end

UTU=0.25*(UTU+UTU');

end