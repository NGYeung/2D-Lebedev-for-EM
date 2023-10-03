function [Out] = InverseBlockMatrix(IN,N1,N2)
%INVERSEBLOCKMATRIX computes the inverse of a 2x2 block diagonal matrix
% off diagonal elements are allowed to be prolongation matrices as long as
% C'*C is diagonal
% In - input of diagonal blokc matrix [A, C, C', B]
% N1 size of A
% N2 size of B

if(length(IN)~= N1+N2)
    error('Incompatible input N1 + N2 not N')
end

%extract principal matricees
A=IN(1:N1,1:N1);
B=IN(N1+(1:N2),N1+(1:N2));
C=IN((1:N1),N1+(1:N2));

A_d=diag(A);
B_d=diag(B);
% C_d_square=diag(C'*C);

Ainv=spdiags(1./A_d,0,N1,N1);
Binv=spdiags(1./B_d,0,N2,N2);

Maa=spdiags(1./(A_d-diag(C*Binv*C')),0,N1,N1);
Mbb=spdiags(1./(B_d-diag(C'*Ainv*C)),0,N2,N2);


Out=[Maa -Maa*C*Binv; -Mbb*C'*Ainv Mbb];

end

