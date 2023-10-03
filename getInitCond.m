function [E_Init,Q,lambda0ftr,StartBox,EndBox,sizes]=getInitCond(src,sizesIn)
%function E_Init=getInitCond(src,sizes)

N=30; %we take 20 point round the source locations
disp('Source function assumes homoegenious vaccum background right now')

E_Init=zeros(sizesIn.NE,src.nsrc);
A_prev=[];
for it=1:src.nsrc
    
    
    %wall detection
    StartBox=max(src.src_loc(1:2,it)-N,0);
    EndBox=min(src.src_loc(1:2,it)+N,[sizesIn.Ny, sizesIn.Nx]');
    
    relativeLoc=min(src.src_loc(1:2,it)-N,0)+N;
    
    NxBox=EndBox(2)-StartBox(2);
	NyBox=EndBox(1)-StartBox(1);
    %set medum in background (in the)
    
    [sizes, ~,~, eps, ~]=get_medium('Homo',NxBox,NyBox);
    [DEH_sym,DHE_sym,~]=getSymmetricOperators(sizes,eps);
    
    if(src.src_loc(3,it)==1)%xsource
        Ex_pd=zeros(sizes.ex1);
        Ex_pd(relativeLoc(2),relativeLoc(1))=1;
        Ex_dp=conv2(conv2(Ex_pd,[1/2 1/2]),[1/2;1/2]);
        
        RHS=[Ex_pd(:);zeros(sizes.ey2p,1);Ex_dp(:);zeros(sizes.ey1p,1)];
    elseif(src.src_loc(3,it)==2)%ysource
        Ey_pd=zeros(sizes.ey2);
        Ey_pd(relativeLoc(2),relativeLoc(1)+1)=1;
        Ey_dp=conv2(conv2(Ey_pd,[1/2 1/2],'valid'),[1/2;1/2],'valid');
        
        RHS=[zeros(sizes.ex1p,1);Ey_pd(:);zeros(sizes.ex2p,1);Ey_dp(:)];
    else
        error('specified source type not supported')
    end
    
    %Computing the source via second order operator
    AA=DHE_sym*DEH_sym;
    AA=-0.5*(AA+AA');

    if(sum(size(A_prev)==size(AA))==2 ) 
        if(norm(A_prev-AA,'fro')/norm(AA,'fro') >1e-10)
            %operators are nto the saem, recompute
            tic;[Q,lambda0ftr]=eig(full(AA),'vector');toc
        else
            disp('Eigendecomposition reused when computing source')
        end
    else
        % operators have different size, recompute
        tic;[Q,lambda0ftr]=eig(full(AA),'vector');toc
    end

%     tic;[Q,E]=eig(full(A));toc % 263 secs
    bSmall=Q * spdiags(sqrt(src.F(sqrt(lambda0ftr))), 0, sizes.NE, sizes.NE) *(Q'*RHS);
                    
%     b=Q*diag(sqrt(src.F(diag(E))))*(Q\RHS);
    
    %embedd back
    mask_ex1=zeros(sizesIn.ex1);
    mask_ey2=zeros(sizesIn.ey2);
    mask_ex2=zeros(sizesIn.ex2);
    mask_ey1=zeros(sizesIn.ey1);

    mask_ex1(StartBox(2)+(1:sizes.ex1(1)),StartBox(1)+(1:sizes.ex1(2)))=reshape(bSmall(sizes.index_ex1p-sizes.NH),sizes.ex1);
    mask_ey2(StartBox(2)+(1:sizes.ey2(1)),StartBox(1)+(1:sizes.ey2(2)))=reshape(bSmall(sizes.index_ey2p-sizes.NH),sizes.ey2);
    mask_ex2(StartBox(2)+(1:sizes.ex2(1)),StartBox(1)+(1:sizes.ex2(2)))=reshape(bSmall(sizes.index_ex2p-sizes.NH),sizes.ex2);
    mask_ey1(StartBox(2)+(1:sizes.ey1(1)),StartBox(1)+(1:sizes.ey1(2)))=reshape(bSmall(sizes.index_ey1p-sizes.NH),sizes.ey1);
    
   
    E_Init(:,it)=[mask_ex1(:);mask_ey2(:);mask_ex2(:);mask_ey1(:)];
    A_prev=AA;
end

    E_Init = 20.*E_Init;

end