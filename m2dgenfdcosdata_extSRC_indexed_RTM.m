function [ imgG1x, imgG1y, imgG1] = m2dgenfdcosdata_extSRC_indexed_RTM( t, A, B,gcosSC_dt,U0,indexV,sizes,show) 


%M2DGENCOSDATA Generate the data for the cosine fitting problem
% from the operator matrix A, source/receiver matrix B  %% WHAT IS
% SOURCE-RECEIVER MATRIX
% using time stepping for both the matrix exponential in the 
% sources/receivers and the cosine itself
% Alexander Mamonov, University of Houston, 2013, 2014, 2015
% Jorn Zimmerling, UMICH, added wavelet support for alex 2019
%==========================================================================
% In this block version the wavelet has the dimensions
% sources x nr of experiments x times
% optimize by computing the wavelet differences once

%

%tau=timingST.tau;
%ntau=timingST.ntau; %#smaller time steps per tau
%t=timingST.ht; %tau (why the rep
%nt=timingST.nt; %number of total timestep
nts = 20;
tau = t/nts;
m1=1.*round(0.9*sizes.Nx/t); %%%%%With or without *2 %0112
%nts=timingST.nts; %smaller time steps per tau
%m1=src.nsrc; % number of source x polarization


n=size(B,2);

ind = @(j) n*(j-1) + (1:n); %This is for find the block in U

% compute the matrix exponential exp(0.5 * tau * A) * B
% as the solution u(0.5 * tau) of a parabolic equation u_t = Au, u(0) = B

% expAtB = B;
% htau = 0.5 * tau / ntau;
%for j = 1:ntau
%    expAtB = expAtB + htau * (A * expAtB);
%end

%==========================================================================
% compute the solution and the data
%%%%Makesure indexV is for the alpha grid only.


hts2 = tau^2;


%----------DONE---------------YIYANG

nr_tr=size(gcosSC_dt,2);%number of parallel experiments  %number of sources in our case?
lengW=size(gcosSC_dt,3); % the number of U0-1 because of the difference

% correct initial conditions
% compatible with the cosine data model
% VERY IMPORTANT TO GET THIS PART RIGHT!
Uprev = zeros(size(B,1),nr_tr); %for now, later more

% one application of the DISCRETE PROPAGATOR
Ucurr = Uprev + (0.5 * hts2) * (A * Uprev);


%% Need two different images.
imgG1x=zeros(size(indexV.x,1),1);
imgG1y=zeros(size(indexV.y,1),1);
% imgG2=zeros(size(XrangeIndex));
% imgG3=zeros(size(XrangeIndex));

iter=1;

% gcos=zeros(size(B,2),nr_tr,nts*nt,'gpuArray'); %??
%minus two as there is one step overlap in recording
for k = 0:((2*m1-1)*nts) %Why 2*m1-1 !!!!! MODIFIED.let's see
%     if mod(k, nts) == 0
%        gcos{(k/nts) + 1}  = B' * Uprev;
%     end
    
    % Iw is the difference of the wavelet integrated over the current and
    % previous timestep
    % I define my wavelet on the dual timesteps and therefore this integral
    % by midpoint rule is hts*(wavelet(k)-wavelet(k-1))
    % If I define the wavelet on the primary time sets 
%     if(k+2<=lengW) 
%        Iw=hts*(wavelet(:,k+2)-wavelet(:,k+1)); % sort of the same...

    if(k+1<=lengW) %changed



        Iw=tau*(gcosSC_dt(:,:,k+1)); %changed 0109
        % Iw=hts2.*(gcosSC_dt(:,:,k+1)); 
        %Iw=hts*(wavelet(:,k+1));
        
        Unext = 2*Ucurr - Uprev + hts2 * (A * Ucurr)+B*Iw; %.*??  %changed + to - 0119
       % Unext = 2*Ucurr - Uprev + hts2 * (A * Ucurr)+expAtB*Iw;  A*Fe_curr+2.*Fe_curr - Fe_prev
    else
        Unext = 2*Ucurr - Uprev + hts2 * (A * Ucurr);
    end
    % If I define the wavelet on the primary time sets I get 
%     Iw=hts*(wavelet(k+1)-wavelet(k-1))
    if(mod(k,nts)==0)
        %integer timestep
        indexNr=(k/nts)
        %xrange index versus r index
        %g is source index , backwards in time, randeindex

        %% Rewrite the index for the Lebedev Grid. Use the alpha grid.
         imgG1x=imgG1x+sum(Uprev(indexV.x,:).*U0(indexV.x,ind(2*m1-indexNr)),2); %%%0105 deleted the *2 change the indexing of U0, =1
         imgG1y=imgG1y+sum(Uprev(indexV.y,:).*U0(indexV.y,ind(2*m1-indexNr)),2);
         imgG1 = imgG1x+ imgG1y;
if show
 figure(11)
subplot(211)
        kp=reshape(Uprev(1:sizes.ex1p,5),sizes.ex1(1),[]);
        imagesc(kp)
        colorbar

        subplot(212)
        kk = B;
        imagesc(reshape(kk(1:sizes.ex1p,5),sizes.ex1(1),[]))
        colorbar
        
    
  %      %axis([1 300 1 150])
%         axis([150])
        
        %set(gcf,'Position',[1000         431        1421         907])
        pause(0.5);    

end        
       
    end
    
    Uprev = Ucurr;
    Ucurr = Unext;
    

  
    
    
end
end
