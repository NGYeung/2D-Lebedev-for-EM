function [Malpha,Mbeta,Malpha_sq,Mbeta_sq,MeInv,Mhinv,Malpha_sq_inv,Mbeta_sq_inv,Me_sq_inv]=MediumMatrix(sizes,eps,mu)
% this functions takes:
% Input:
%       eps (Nx x Ny) and
%       mu  (Nx x Ny)
%
% averages the media to the respective grids and then constructe the mediu
% matrices
% Jorn Zimmerling @ umich 2022

if(any(mu(:) ~= ones(size(mu(:)))))
    warning('no mu inmpe=lemented as of yet')
end

%we assume there is no constract in mu
% MhInv=speye(length(sizes.NH));
Mhinv=speye(length(sizes.NH),length(sizes.NH));

%we average the medium to ex (p x d) and ey (d x p)
% we assume that the outmost roaw and columns f eps are a homogenious meidum
for it=1:3
    eps_pd_til(:,:,it)=conv2(squeeze([ eps(:,1,it) eps(:,:,it) eps(:,end,it)]),[1/2 1/2],'valid');
    eps_dp_til(:,:,it)=conv2(squeeze([ eps(1,:,it);eps(:,:,it);eps(end,:,it)]),[1/2;1/2],'valid');
    
    eps_pd_til_tr(:,:,it)=eps_pd_til(:,2:(end-1),it);
    eps_dp_til_tr(:,:,it)=eps_dp_til(2:(end-1),:,it);
end


%Compute squareroot of meedium matrix - tested and checked
eps_pd_til_insq   =InverseMat(SquareRootMat(eps_pd_til));
eps_dp_til_insq   =InverseMat(SquareRootMat(eps_dp_til));
eps_pd_til_tr_insq=InverseMat(SquareRootMat(eps_pd_til_tr));
eps_dp_til_tr_insq=InverseMat(SquareRootMat(eps_dp_til_tr));

eps_pd_til_sq   =(SquareRootMat(eps_pd_til));
eps_dp_til_sq   =(SquareRootMat(eps_dp_til));
eps_pd_til_tr_sq=(SquareRootMat(eps_pd_til_tr));
eps_dp_til_tr_sq=(SquareRootMat(eps_dp_til_tr));


% alpha grid primary x dual  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=eps_pd_til_tr(:,:,1);
Mepsx=spdiags(ee(:),0,sizes.ex1p,sizes.ex1p);
Mepsxinv=spdiags(1./ee(:),0,sizes.ex1p,sizes.ex1p);

ee=eps_pd_til(:,:,2);
Mepsy=spdiags(ee(:),0,sizes.ey2p,sizes.ey2p);
Mepsyinv=spdiags(1./ee(:),0,sizes.ey2p,sizes.ey2p);

ee=eps_pd_til(:,:,3);
Mepsxy=spdiags(ee(:),0,sizes.ey2p,sizes.ey2p);
mask=logical([zeros(size(eps(:,1,1))) ones(size(eps_pd_til_tr(:,:,1))) zeros(size(eps(:,1,1)))] ) ;
Mepsxy=Mepsxy(mask(:),:);


Malpha=[Mepsx Mepsxy;Mepsxy' Mepsy];

Maa=spdiags(1./diag(Mepsx-Mepsxy*Mepsyinv*Mepsxy'),0,sizes.ex1p,sizes.ex1p);% precompute this on the right grid first?
Mbb=spdiags(1./diag(Mepsy-Mepsxy'*Mepsxinv*Mepsxy),0,sizes.ey2p,sizes.ey2p);

MaInv=...
[Maa -Maa*Mepsxy*Mepsyinv;...
 -Mbb*Mepsxy'*Mepsxinv Mbb];

% beta grid dual x primary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=eps_dp_til(:,:,1);
Mepsx=spdiags(ee(:),0,sizes.ex2p,sizes.ex2p);
Mepsxinv=spdiags(1./ee(:),0,sizes.ex2p,sizes.ex2p);

ee=eps_dp_til_tr(:,:,2);
Mepsy=spdiags(ee(:),0,sizes.ey1p,sizes.ey1p);
Mepsyinv=spdiags(1./ee(:),0,sizes.ey1p,sizes.ey1p);

%mask matricees
%delete zeros at boundary for oppdiagonal ones
ee=eps_dp_til(:,:,3);
Mepsxy=spdiags(ee(:),0,sizes.ex2p,sizes.ex2p);
mask=logical([zeros(size(eps(1,:,1))); ones(size(eps_dp_til_tr(:,:,1))); zeros(size(eps(1,:,1)))] );
Mepsxy=Mepsxy(mask(:),:)';

Mbeta=[Mepsx Mepsxy;Mepsxy' Mepsy];

Maa=spdiags(1./diag(Mepsx-Mepsxy*Mepsyinv*Mepsxy'),0,sizes.ex2p,sizes.ex2p);% precompute this on the right grid first?
Mbb=spdiags(1./diag(Mepsy-Mepsxy'*Mepsxinv*Mepsxy),0,sizes.ey1p,sizes.ey1p);

MbInv=...
[Maa -Maa*Mepsxy*Mepsyinv;...
 -Mbb*Mepsxy'*Mepsxinv Mbb];

MeInv=...
[MaInv  sparse(sizes.Nalpha,sizes.Nbeta);...
 sparse(sizes.Nalpha,sizes.Nbeta)  MbInv ];



% inverseSquareRoots

% alpha grid primary x dual  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=eps_pd_til_tr_insq (:,:,1);
Mepsx_insq=spdiags(ee(:),0,sizes.ex1p,sizes.ex1p);

ee=eps_pd_til_insq (:,:,2);
Mepsy_insq =spdiags(ee(:),0,sizes.ey2p,sizes.ey2p);

ee=eps_pd_til_insq (:,:,3);
Mepsxy_insq =spdiags(ee(:),0,sizes.ey2p,sizes.ey2p);
mask=logical([zeros(size(eps(:,1,1))) ones(size(eps_pd_til_tr(:,:,1))) zeros(size(eps(:,1,1)))] ) ;
Mepsxy_insq=Mepsxy_insq(mask(:),:);

Malpha_sq_inv=[Mepsx_insq Mepsxy_insq;Mepsxy_insq' Mepsy_insq];

% beta grid dual x primary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=eps_dp_til_insq(:,:,1);
Mepsx_insq=spdiags(ee(:),0,sizes.ex2p,sizes.ex2p);

ee=eps_dp_til_tr_insq(:,:,2);
Mepsy_insq=spdiags(ee(:),0,sizes.ey1p,sizes.ey1p);

%mask matricees
%delete zeros at boundary for oppdiagonal ones
ee=eps_dp_til_insq(:,:,3);
Mepsxy_insq=spdiags(ee(:),0,sizes.ex2p,sizes.ex2p);
mask=logical([zeros(size(eps(1,:,1))); ones(size(eps_dp_til_tr(:,:,1))); zeros(size(eps(1,:,1)))] );
Mepsxy_insq=Mepsxy_insq(mask(:),:)';



Mbeta_sq_inv=[Mepsx_insq Mepsxy_insq;Mepsxy_insq' Mepsy_insq];

Me_sq_inv=...
[Malpha_sq_inv sparse(sizes.Nalpha,sizes.Nbeta);...
sparse(sizes.Nbeta,sizes.Nalpha) Mbeta_sq_inv];
%%%%%%%%%%%%%%%% square roots
% alpha grid primary x dual  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=eps_pd_til_tr_sq (:,:,1);
Mepsx_sq=spdiags(ee(:),0,sizes.ex1p,sizes.ex1p);

ee=eps_pd_til_sq (:,:,2);
Mepsy_sq =spdiags(ee(:),0,sizes.ey2p,sizes.ey2p);

ee=eps_pd_til_sq (:,:,3);
Mepsxy_sq =spdiags(ee(:),0,sizes.ey2p,sizes.ey2p);
mask=logical([zeros(size(eps(:,1,1))) ones(size(eps_pd_til_tr(:,:,1))) zeros(size(eps(:,1,1)))] ) ;
Mepsxy_sq=Mepsxy_sq(mask(:),:);

Malpha_sq=[Mepsx_sq Mepsxy_sq;Mepsxy_sq' Mepsy_sq];

% beta grid dual x primary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee=eps_dp_til_sq(:,:,1);
Mepsx_sq=spdiags(ee(:),0,sizes.ex2p,sizes.ex2p);

ee=eps_dp_til_tr_sq(:,:,2);
Mepsy_sq=spdiags(ee(:),0,sizes.ey1p,sizes.ey1p);

%mask matricees
%delete zeros at boundary for oppdiagonal ones
ee=eps_dp_til_sq(:,:,3);
Mepsxy_sq=spdiags(ee(:),0,sizes.ex2p,sizes.ex2p);
mask=logical([zeros(size(eps(1,:,1))); ones(size(eps_dp_til_tr(:,:,1))); zeros(size(eps(1,:,1)))] );
Mepsxy_sq=Mepsxy_sq(mask(:),:)';

Mbeta_sq=[Mepsx_sq Mepsxy_sq;Mepsxy_sq' Mepsy_sq];

end

% auxiliarry functions to manipulate our summetric tensor defition
function eps_sq=SquareRootMat(eps)
% square root or positive definite symmetric tensor
%
%[1 3;
% 3 2]
%

eps_sq=eps;
s=sqrt(eps(:,:,1).*eps(:,:,2)-eps(:,:,3).^2);
t=sqrt(eps(:,:,1)+eps(:,:,2)+2*s);
eps_sq(:,:,1)=(eps_sq(:,:,1)+s)./t;
eps_sq(:,:,2)=(eps_sq(:,:,2)+s)./t;
eps_sq(:,:,3)=eps_sq(:,:,3)./t;

end

function eps_inv=InverseMat(eps)

eps_sq=eps;
s=eps(:,:,1).*eps(:,:,2)-eps(:,:,3).^2;

eps_inv(:,:,1)=(eps_sq(:,:,2))./s;
eps_inv(:,:,2)=(eps_sq(:,:,1))./s;
eps_inv(:,:,3)=-eps_sq(:,:,3)./s;

end

function eps_mult=MultMat(eps)

eps_mult=zeros(size(eps));
eps_mult(:,:,1)=eps(:,:,1).^2+eps(:,:,3).^2;
eps_mult(:,:,2)=eps(:,:,2).^2+eps(:,:,3).^2;
eps_mult(:,:,3)=(eps(:,:,2)+eps(:,:,1)).*eps(:,:,3);

end