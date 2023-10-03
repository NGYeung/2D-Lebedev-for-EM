function [U_ave] = e_field_averaging(U,sizes)
%Fix the input tomorrow
%This function average the e-field to hz1 grid points (all dual). 
%This function also constructs the norm



UEX1 = U(sizes.index_ex1p-sizes.NH,:);
UEX2 = U(sizes.index_ex2p-sizes.NH,:);
UEY1 = U(sizes.index_ey1p-sizes.NH,:);
UEY2 = U(sizes.index_ey2p-sizes.NH,:);

tm = size(U,2);




%reshape
UEX1 = reshape(UEX1,sizes.ex1(1),sizes.ex1(2),tm);
UEX2 = reshape(UEX2,sizes.ex2(1),sizes.ex2(2),tm);
UEY1 = reshape(UEY1,sizes.ey1(1),sizes.ey1(2),tm);
UEY2 = reshape(UEY2,sizes.ey2(1),sizes.ey2(2),tm);

U_ave = zeros(2*sizes.hz1p ,tm);

%averaging [each time slice.] Write this shit tomorrow in a loop.


for t = 1:tm

UEXt = 1.*conv2(UEX1(:,:,t),[1/2 1/2])+1.*conv2(UEX2(2:(end-1),:,t),[1/2 1/2]');%modified
UEYt = 1.*conv2(UEY1(:,:,t),[1/2 1/2]')+1.*conv2(UEY2(:,2:(end-1),t),[1/2 1/2]);%
UEXt = reshape(UEXt,sizes.hz1p ,1);
UEYt = reshape(UEYt,sizes.hz1p ,1);

U_ave(1:sizes.hz1p ,t) = UEXt;
U_ave(sizes.hz1p +(1:sizes.hz1p),t) = UEYt;

end



end

%and then we can do the V = U\R outside the function
