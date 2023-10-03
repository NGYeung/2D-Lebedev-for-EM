function [sizes, timing,src, eps, mu]=get_medium(medium,varargin)
% This functions provides some testcases for FDTD
%
% input: medium is a string
%
% output: sizes - structure that contains all size informatoin of the
%                 domain and the fields
%         timing - timing file for FDTD
%         src    - source information
%         eps    - eps Nx x Ny x 3 is the symemtric tensor [exx, eyy, exy]
%         mu     - mu is the scalar magnetic medium

   Nx=varargin{1};
   Ny=varargin{2};
   

if(strcmp(medium,'testAnisotropic'))

   sizes.delta=1;
   
   
   sizes.Nx=Nx;
   sizes.Ny=Ny;
   
   timing=[];
   src=[];
   
   
   mu=ones(sizes.Nx,sizes.Ny);
   eps=ones(sizes.Nx,sizes.Ny,3);
   eps(:,:,3)=0;
%    eps(:,:,3)=1.5; % epsxy
%    eps(:,:,1)=2.5; % exx
%    eps(:,:,2)=2.5; % exy
%    %direction 1 1/1 -1 with eps=1 4
%    
   
   eps(50:60, 50:70,3)=1.5; % epsxy
   eps(60:64,50:70,1)=2.5; % exx
   eps(60:64,50:70,2)=2.5; % exy
   %direction 1 1/1 -1 with eps=1 4
   
   
   
   
   
elseif(strcmp(medium,'testAnisotropicXY'))

   sizes.delta=1;
   
   
   sizes.Nx=Nx;
   sizes.Ny=Ny;
   
   timing=[];
   src=[];
   
   
   mu=ones(sizes.Nx,sizes.Ny);
   eps=ones(sizes.Nx,sizes.Ny,3);
   eps(:,:,3)=0;
%    eps(:,:,3)=1.5; % epsxy
%    eps(:,:,1)=2.5; % exx
%    eps(:,:,2)=2.5; % exy
%    %direction 1 1/1 -1 with eps=1 4
%    
   
   eps(60:64, 50:70,3)=2; % epsxy
   eps(60:64,50:70,1)=1.5; % exx
   eps(60:64,50:70,2)=1.5; % exy
   %direction 1 1/1 -1 with eps=1 4


elseif(strcmp(medium,'testAnisotropicX'))

   sizes.delta=1;
   
   
   sizes.Nx=Nx;
   sizes.Ny=Ny;
   
   timing=[];
   src=[];
   
   
   mu=ones(sizes.Nx,sizes.Ny);
   eps=ones(sizes.Nx,sizes.Ny,3);
   eps(:,:,3)=0;
%    eps(:,:,3)=1.5; % epsxy
%    eps(:,:,1)=2.5; % exx
%    eps(:,:,2)=2.5; % exy
%    %direction 1 1/1 -1 with eps=1 4
%    
   
  % eps(60:64, 50:70,3)=2; % epsxy
   eps(40:64,50:70,1)=8; % exx
   eps(40:64,50:70,2)=1; % exy
   %direction 1 1/1 -1 with eps=1 4


elseif(strcmp(medium,'testAnisotropicY'))

   sizes.delta=1;
   
   
   sizes.Nx=Nx;
   sizes.Ny=Ny;
   
   timing=[];
   src=[];
   
   
   mu=ones(sizes.Nx,sizes.Ny);
   eps=ones(sizes.Nx,sizes.Ny,3);
   eps(:,:,3)=0;
%    eps(:,:,3)=1.5; % epsxy
%    eps(:,:,1)=2.5; % exx
%    eps(:,:,2)=2.5; % exy
%    %direction 1 1/1 -1 with eps=1 4
%    
   
  % eps(60:64, 50:70,3)=2; % epsxy
   eps(60:64,50:70,1)=1; % exx
   eps(60:64,50:70,2)=2; % exy
   %direction 1 1/1 -1 with eps=1 4


elseif(strcmp(medium,'Homo'))
   Nx=varargin{1};
   Ny=varargin{2};
   
   sizes.delta=1;
   
   
   sizes.Nx=Nx;
   sizes.Ny=Ny;
   
   timing=[];
   src=[];
   
   mu=ones(sizes.Nx,sizes.Ny);
   eps=ones(sizes.Nx,sizes.Ny,3);
   eps(:,:,3)=0;

elseif(strcmp(medium(1:6),'config'))
  Str = strcat('medium_lib/',medium) ;
  run(Str)


else 
    error(['Medium ' medium ' not defined in get_medium.m '])
end




% sizes in each cluster
sizes.hz1 = [Nx   Ny];
sizes.hz2 = [Nx+1   Ny+1]; % exclude outer boundary for Neumann!!
sizes.ex1 = [Nx   Ny-1];
sizes.ey2 = [Nx   Ny+1];  % if ey is dirichlet!
sizes.ey1 = [Nx-1 Ny];
sizes.ex2 = [Nx+1 Ny];  % if ex is dirichlet!

sizes.hz1p = prod(sizes.hz1);
sizes.hz2p = prod(sizes.hz2);
sizes.ex1p = prod(sizes.ex1);
sizes.ey2p = prod(sizes.ey2);
sizes.ey1p = prod(sizes.ey1);
sizes.ex2p = prod(sizes.ex2);

sizes.index_hz1p = 1:sizes.hz1p;
sizes.index_hz2p = sizes.index_hz1p(end)+(1:sizes.hz2p);
sizes.index_ex1p = sizes.index_hz2p(end)+(1:sizes.ex1p);
sizes.index_ey2p = sizes.index_ex1p(end)+ (1:sizes.ey2p);
sizes.index_ex2p = sizes.index_ey2p(end)+(1:sizes.ex2p);
sizes.index_ey1p = sizes.index_ex2p(end)+(1:sizes.ey1p);


sizes.NH     =sizes.hz1p+sizes.hz2p;
sizes.Nalpha =sizes.ex1p+sizes.ey2p;
sizes.Nbeta  =sizes.ey1p+sizes.ex2p;
sizes.NE     =sizes.Nalpha+sizes.Nbeta;

sizes.Ntot=sizes.NH + sizes.Nalpha + sizes.Nbeta;

%getGridForPlotting - delta=1 - starting at 0
%truncated grid sizes
sizes.grids.ex1{1}=0.5:1:(Nx-0.5);sizes.grids.ex1{2}=1:(Ny-1);
sizes.grids.ey1{1}=1:1:(Nx-1);sizes.grids.ey1{2}=0.5:(Ny-0.5);
sizes.grids.hz1{1}=0.5:1:(Nx-0.5);sizes.grids.hz1{2}=0.5:(Ny-0.5);

sizes.grids.ex2{1}=0:Nx;sizes.grids.ex2{2}=0.5:(Ny-0.5);
sizes.grids.ey2{1}=0.5:1:(Nx-0.5);sizes.grids.ey2{2}=0:Ny;
sizes.grids.hz2{1}=0:(Nx);sizes.grids.hz2{2}=0:Ny;
end