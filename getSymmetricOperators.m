function [DEH_sym,DHE_sym,A]=getSymmetricOperators(sizes,eps)
%function [DEH_sym,DHE_sym,A]=getSymmetricOperators(sizes,eps,mu)
%
% Constructs first order operator A (antisymetric)and computents DEH_sym
% and DHE_sym
% 
%

%disp('Homogenious mu expected')

[Me_sq_inv]=getMediumMatrixEPS(sizes,eps);

%operators
[~,~,~,~,DEH,DHE,~,~,WH_sq,WE_sq,WH_sq_inv,WE_sq_inv]=get_D_2D(sizes);



% writeup the magical symmetries, prolongation operators and in to doc
%% Construct symmetrized operator
DEH_sym=WH_sq*DEH*WE_sq_inv*Me_sq_inv;
DHE_sym=Me_sq_inv*WE_sq*DHE*WH_sq_inv;
%symmetry - checked
% norm(DEH_sym.'+DHE_sym,'F')
A=[sparse(sizes.NH,sizes.NH) DEH_sym;
   DHE_sym sparse(sizes.NE,sizes.NE)];
A=0.5*(A-A.');%symmetrize



end