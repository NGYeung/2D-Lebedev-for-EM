
function [] = RUN_Lebedev(RTMflag, Nx, Ny, medium_id, folder, sourceflag, tau, sep, thetaset)
% Input 1: boolean, RTMflag = 0 --> NormG    RTMflag = 1 --> RTM
% Input 2, 3: integers, Nx, Ny, the sizes of the finite difference grid
% Input 4: string, medium_id, link to the medium in mediumlib
% Input 5: string, The folder in which results are stored
% Input 6: boolean, sourceflag = 0 use presaved source       sourflag = 1 compute new source
% Input 7: double, timestep for samples
% Input 8: integer, seperation between sources
% Input 9: 2x2 matrix, each column being a unit vector that indicates the direction (for imaging)




fname = folder;

%imagerange.x=40:Nx*2-40;%without averaging
%imagerange.y=40:Ny*2-40;%without averaging


imagerange.x=20:Nx-20;%without averaging %change it back to 100 pls %%%%0105 remember to switch back
imagerange.y=20:Ny-20;%without averaging

%imagerange.x=20:Nx-20;%without averaging %change it back to 100 pls
%imagerange.y=20:Ny-20;%without averaging

set(0,'defaulttextInterpreter','latex') %latex axis labels
set(0, 'defaultAxesTickLabelInterpreter','latex'); set(0, 'defaultLegendInterpreter','latex');

%%
medium = medium_id;
[sizes, timing,src, eps, mu]=get_medium(medium, Nx,Ny, sourceflag);
medium0='Homo';[sizes, timing,src, eps0, mu0]=get_medium(medium0,Nx,Ny);


[DEH_sym,DHE_sym,A]=getSymmetricOperators(sizes,eps);  %%%%construct A the right way.
[DEH_sym0,DHE_sym0,A0]=getSymmetricOperators(sizes,eps0);

%compute initial condition - source information
% 15 ==> 75 cutoff: 8 => /40; cutoff: 12 ==>/60  cutoff: 6 => 30
src.wc = 6*pi/60; 
src.wb = 4*pi/60;
src.wcut=src.wc+src.wb;
src.c0=1;src.delta=1;   %should be routed from code and not set hard
src.lambda_cut=src.c0*(2*pi/src.wcut);
src.ppw_cut=src.lambda_cut/src.delta;

src.F = @(omega) omega.^2.*(0.5*( exp(- ((omega - src.wc).^2)./((2*src.wb ).^2) ) ...
                 + exp(- ((omega + src.wc).^2)./((2*src.wb ).^2) ) ));

%%**
src.nsrc=2*round(Ny-40)/sep; %changed to 6 instead of 4 at 11/18


% SOURCE LOCATION
src.src_loc=[round(linspace(20,sizes.Ny-20,src.nsrc/2)) round(linspace(20,sizes.Ny-20,src.nsrc/2)); 4*ones(1,src.nsrc);  %changed 0112
             ones(1,src.nsrc/2) 2*ones(1,src.nsrc/2)];








%%
% calculate wavelength for plotting use.
% central wavelength
disp('wavelength at cut-off: ')
lamb_c = 2*pi/src.wcut
%express tau in terms of wc:
% T = 2pi/src.wcut
% Nyquist rate = pi/src.wcut
Nyqst = pi/src.wcut;
r = tau/Nyqst;
lamb_c = 2*pi/src.wc;

%%
         
filename = strcat('sources/E_Init_',num2str(Nx*Ny),'_',num2str(src.nsrc),'.mat');       
         
         
if sourceflag == true

[E_Init,Q,Eig]=getInitCond(src,sizes);

save(filename, 'E_Init');

end

if sourceflag == false

E_Init = importdata(filename);
%Q = importdata('Q.mat');
%Eig = importdata('eigenvalue.mat');
end

%%

%make sure initial condition in real
if(isreal(E_Init(:)) ==0 )
    warning('Initial condition not real. Is wavelet even or rounding error?')
   if(norm(real(E_Init(:)))/norm(E_Init(:))>0.99)
       disp('Initial condition mostly real')
       E_Init=real(E_Init);
   else
       error('Check initial conditions')
   end
    
end


%%****
run('Config')


%%
show= 0;
disp('Run Forward model on CPU')
tic;



%%
%------------------------------------RTM --------------------------------------------------------------------------------------
if RTMflag

[D,U_E,~,RTMD]=timestep_sym_init_2nd_test(A,E_Init,show,sizes, tau, 3,imagerange,src); %second order

toc;


show = 0;
disp('Run Forward model - background')
tic;


[D0,U_E0,Ac0,RTMD0]=timestep_sym_init_2nd_test(A0,E_Init,show,sizes, tau, 3, imagerange,src);


toc;
%%

D_sc = RTMD-RTMD0; 
D_sc = D_sc(:,:,size(D_sc,3):-1:1);
D_sc_dt = diff(D_sc,1,3);
clear U_sc

%---------------------------------------------------------------------------------------------------------------------------------


else


run('Config')

%%
%----------------------------------------------------------NORMG FORWARD-----------------------------------------------------------

[D,U_E,Ac]=timestep_sym_init_2nd_test(A,E_Init,show,sizes, tau, 1,imagerange,src);

toc;


%%

show = 0;
disp('Run Forward model - background')
tic;

[D0,U_E0,Ac0]=timestep_sym_init_2nd_test(A0,E_Init,show,sizes, tau, 1, imagerange,src);
toc;


%%
sizes

%Plot Data

%PlotData



end


%%
%-----------------------------------------------------IMAGING SCRIPT-----------------------------------------------------------------------------


run('Run_Image')




end