function [Data]=timestep_sym_init_2nd_gpu(A,E_Init,show,sizes,t,cluster,imagerange,src)



cluster = 1;


nt=20;
tau=t/nt;
%nts=round(max(min(1.8*t,1.4),1.2)*sizes.Nx/tau/nt)
nts=2*round(0.7*sizes.Nx/tau/nt); %0109 remember to change it
% to make sure that the wave from the bottom doesn't reach the receiver but penetrate deep enough 

MAXStep=nts*nt;
A=(tau^2).*A(sizes.NH+(1:sizes.NE),:)*A(:,sizes.NH+(1:sizes.NE));
A = gpuArray(A);
%disp('update!')
%tic;
Fe_curr=E_Init;
Fe_prev = 0.5*A*Fe_curr+Fe_curr; %%%%%check on the update rules!!!!!
% Fe_prev = A*Fe_curr;
%toc;

%Fh_curr=+0.5*tau*D*Fe_curr;
% Fh_next=zeros(size(Fh_curr));
% Fe_next=zeros(size(Fe_curr));

nsrc=size(E_Init,2);

%output and storing matrices
Data = zeros(nsrc,nsrc,nts, 'gpuArray');

RTMData = zeros(nsrc,nsrc,MAXStep, 'gpuArray');  %Remember to switch back
%-----------------------------------------------------


%----------------------------------------------------- 

%blockindexing
indBl = @(x)  (1:nsrc)+(x-1)*nsrc;
storID=1;

%----------------------------------------------------- 

if cluster == 1
    IndexCurr=[sizes.index_ex1p sizes.index_ey2p]-sizes.NH;

elseif cluster == 0
    xx = size(imagerange.x,2);
    yy = size(imagerange.y,2);

    maskx = zeros(sizes.ex1(1),sizes.ex1(2));
    maskx(imagerange.x,imagerange.y)=1;
    indexx = find(maskx(:));

    masky = zeros(sizes.ey2(1),sizes.ey2(2));
    masky(imagerange.x,imagerange.y)=1;
    indexy = find(masky(:))+sizes.ex1p;

    clear maskx
    clear masky

    IndexCurr=[indexx indexy];
    
    else
    IndexCurr=1:size(E_Init,1);
end
%Data(:,:,storID)=E_Init(IndexCurr,:)'*Fe_curr(IndexCurr,:);

script_deltafun

Data(:,:,storID) = E_Init'*Fe_curr;

%----------------------------------------------------- 


for it=2:MAXStep %

%disp('update!')
%tic
    Fe_next=A*Fe_curr+2.*Fe_curr - Fe_prev; %step k

    Fe_prev = Fe_curr;
    Fe_curr=Fe_next;
%toc
 
    % compute data
    if floor((it-1)/nt)==(it-1)/nt   %changed the old one should be in the backup file.
%     disp('Save Data!')
  
        storID=(it-1)/nt+1;
        %storID=it/nt+1;
        
        Data(:,:,storID)=E_Init'*Fe_curr;
    end
    
    
    if(and(show,floor(it/10)==it/10))
        %probably I shoudl construct the grid to indicate the shifts where
        %thy live
        Fst=[zeros(sizes.NH,nsrc);Fe_curr];
        figure(6)
         set(gcf,'Position',[100         100        875         815])
        %plot hz
            subplot(231)
                PlField=reshape(Fst(sizes.index_hz1p),sizes.hz1);
                imagesc(PlField)
                plotSet(PlField)
                title('$H_z (00)$')
            subplot(234)
                PlField=reshape(Fst(sizes.index_hz2p),sizes.hz2);
                imagesc(PlField)
                plotSet(PlField)
                title('$H_z (11)$')
        %plot ex
            subplot(232)
                PlField=reshape(Fst(sizes.index_ex1p),sizes.ex1);
                imagesc(PlField)
                plotSet(PlField)
                title('$E_x (00)$')
            subplot(235)
                PlField=reshape(Fst(sizes.index_ex2p),sizes.ex2);
                imagesc(PlField)
                plotSet(PlField)
                title('$E_x (11)$')
        %plot rx
            subplot(233)
                PlField=reshape(Fst(sizes.index_ey1p),sizes.ey1);
                imagesc(PlField)
                plotSet(PlField)
                title('$E_y (00)$')
            subplot(236)
                PlField=reshape(Fst(sizes.index_ey2p),sizes.ey2);
                imagesc(PlField)
                plotSet(PlField)
                title('$E_y (11)$')
                
                pause(0.02)
    end
    
end


end

function plotSet(PlField)
crameri('broc'),%caxis([-1 1]*max(abs(PlField(:))));
axis equal, axis tight
colorbar
end