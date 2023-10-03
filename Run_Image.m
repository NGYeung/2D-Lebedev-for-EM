
%%

xx = size(imagerange.x,2);
yy = size(imagerange.y,2);

maskx = zeros(sizes.ex1(1),sizes.ex1(2));
maskx(imagerange.x,imagerange.y)=1;
indexx = find(maskx(:));

masky = zeros(sizes.ey2(1),sizes.ey2(2));
masky(imagerange.x,imagerange.y+1)=1; %changed0314
indexy = find(masky(:));


%%
if RTMflag

indexV.x = indexx; indexV.y = indexy;
%Run the RTM script
%E_Init = E_Init(1:sizes.Nalpha,:);

script_deltafun
[imgG1x, imgG1y, imgG1] = m2dgenfdcosdata_extSRC_indexed_RTM( tau, Ac0, E_Init ,D_sc_dt,U_E0, indexV, sizes,0);  % A is somehow %changeed Dfn
%%
imgG1x = real(reshape(imgG1x, xx,yy));
imgG1y = real(reshape(imgG1y, xx,yy));
imgG1 = real(reshape(imgG1, xx,yy));
%%lam
%%
srtm = figure()
 



k = diff(imgG1);
%k=imgaussfilt(k,[0.5 1e-12]);

 imagesc(k)
 
 axis equal;axis tight;
  %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
ax = gca;
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
      caxis([-1 1]*(max(abs(k(:)))))
       %caxis([-0.1 0.1])
  crameri('broc')
 title('RTM','FontSize',14)
 colorbar



imagetitlefig = strcat('results/', fname, '/', medium_id, '_RTM_' , num2str(tau), '.fig');
imagetitlepng = strcat('results/', fname, '/', medium_id, '_RTM_' , num2str(tau), '.png');
saveas(srtm,imagetitlefig) 
exportgraphics(ax,imagetitlepng,'Resolution',300)
%saveas(srtm,imagetitlepng)
filenamefn = strcat('results/',folder,'/RTM_',medium_id,'_',num2str(sep),'_',num2str(tau));
save(strcat(filenamefn,'imgG1.mat'), 'imgG1');
save(strcat(filenamefn,'imgG1x.mat'), 'imgG1x');
save(strcat(filenamefn,'imgG1y.mat'), 'imgG1y');

%save(strcat(filenamefn,'RTM_IMEx.mat'), 'imgG1x');
%save(strcat(filenamefn,'RTM_IMEx0.mat'), 'imgG1y');





else


PlotTrueSolution


%% Get the orthogonalized basis for both the inclusion and the background.
disp('compute orthogonalized basis')

[R,R0,statOut]=SimStab_withProb_func_tt(D,D0,src,1);


disp('for debug!')
size(R)
size(D)
size(U_E)
src

V = U_E/R;
V0 = U_E0/R0;


%%
%---------------------------Plot Delta Function---------------------------------------------


%PlotDelta_f
%PlotDelta_f

%PlotSolution

xx = size(imagerange.x,2);
yy = size(imagerange.y,2);

maskx = zeros(sizes.ex1(1),sizes.ex1(2));
maskx(imagerange.x,imagerange.y)=1;
indexx = find(maskx(:));

masky = zeros(sizes.ey2(1),sizes.ey2(2));
masky(imagerange.x,imagerange.y+1)=1; %changed0314
indexy = find(masky(:));


%%   
%-------------------------Cleanup------------------------------------

clear U_E
if RTMflag == true
U_dummy = U_E0;
end
clear D

clear D0

ss= src.nsrc/2;
numstep = size(V,2)/src.nsrc;
Uest = V0*R;
Uest = reshape(Uest,size(V,1),src.nsrc,numstep);

Uxx_est = Uest(indexx,1:ss,:);
Uxx_est = reshape(Uxx_est,size(indexx,1),ss*numstep);
Uxy_est = Uest(indexx,ss+(1:ss),:);
Uxy_est = reshape(Uxy_est,size(indexx,1),ss*numstep);
Uyx_est = Uest(sizes.ex1p+indexy,1:ss,:);
Uyx_est = reshape(Uyx_est,size(indexx,1),ss*numstep);
Uyy_est = Uest(sizes.ex1p+indexy,ss+(1:ss),:);
Uyy_est = reshape(Uyy_est,size(indexx,1),ss*numstep);

U_E = V*R;
U_E = reshape(U_E,size(V,1),src.nsrc,numstep);
Uxx_r = U_E(indexx,1:ss,:);
Uxx_r = reshape(Uxx_r,size(indexx,1),ss*numstep);
Uxy_r = U_E(indexx,ss+(1:ss),:);
Uxy_r = reshape(Uxy_r,size(indexx,1),ss*numstep);
Uyx_r = U_E(sizes.ex1p+indexy,1:ss,:);
Uyx_r = reshape(Uyx_r,size(indexx,1),ss*numstep);
Uyy_r = U_E(sizes.ex1p+indexy,ss+(1:ss),:);
Uyy_r = reshape(Uyy_r,size(indexx,1),ss*numstep);

U_E0 = V0*R0;
U_E0 = reshape(U_E0,size(V,1),src.nsrc,numstep);
Uxx_0 = U_E0(indexx,1:ss,:);
Uxx_0 = reshape(Uxx_0,size(indexx,1),ss*numstep);
Uxy_0 = U_E0(indexx,ss+(1:ss),:);
Uxy_0 = reshape(Uxy_0,size(indexx,1),ss*numstep);
Uyx_0 = U_E0(sizes.ex1p+indexy,1:ss,:);
Uyx_0 = reshape(Uyx_0,size(indexx,1),ss*numstep);
Uyy_0 = U_E0(sizes.ex1p+indexy,ss+(1:ss),:);
Uyy_0 = reshape(Uyy_0,size(indexx,1),ss*numstep);



clear U_E
clear U_E0
if RTMflag == true
U_E0 = U_dummy;
end

%%
% 
disp('construct image function')

theta =0;


%length_ = size(imagerange.x,2)*size(imagerange.y,2);
 
 


    Ix1 = sum(Uxx_est.^2,2);
    
    Ixr1 = sum(Uxx_r.^2,2);
    
    Ix01 = sum(Uxx_0.^2,2);


    Iy2 = sum(Uyy_est.^2,2);

    Iyr2 = sum(Uyy_r.^2,2);

    Iy02 = sum(Uyy_0.^2,2);


    Iyx = sum(Uyx_est.^2,2);

    Iyxr = sum(Uyx_r.^2,2);

    Iyx0 = sum(Uyx_0.^2,2);


    Ixy = sum(Uxy_est.^2,2);

    Ixyr = sum(Uxy_r.^2,2);

    Ixy0 = sum(Uxy_0.^2,2);

    







%%


clear maskx
clear masky

IMEA1 = real(reshape(Ix1, xx,yy));
IMEA1r = real(reshape(Ixr1, xx,yy));
IMEA10 = real(reshape(Ix01, xx,yy));
IMEA2 = real(reshape(Iy2, xx,yy));
IMEA2r = real(reshape(Iyr2, xx,yy));
IMEA20 = real(reshape(Iy02, xx,yy));
IMEA3 = real(reshape(Ixy, xx,yy));
IMEA3r = real(reshape(Ixyr, xx,yy));
IMEA30 = real(reshape(Ixy0, xx,yy));
IMEA4 = real(reshape(Iyx, xx,yy));
IMEA4r = real(reshape(Iyxr, xx,yy));
IMEA40 = real(reshape(Iyx0, xx,yy));

%%

s3=figure()
 

k1 = diff(IMEA1-(IMEA10));
%k1 = diff(IMEA1./(IMEA10+1e-14));

 imagesc(k1)
 
 axis equal;axis tight;
  %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
 ax = gca;
 
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 
 
 caxis([-1 1]*max(abs(real(k1(:)))))
 crameri('broc')
 colorbar


Title = strcat('norm g range derivative ($xx$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)


 
imagetitlefig = strcat('results/', fname, '/Ex_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/Ex_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)

 
 
%---------------------------------

 
s3=figure()




 k = diff(IMEA2-(IMEA20));

 imagesc(k)
 
 axis equal;axis tight;
 %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
 ax =gca;
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 caxis([-1 1]*max(abs(real(k(:)))))
  %caxis([-0.1 0.1])
 crameri('broc')
 colorbar


Title = strcat('norm g range derivative ($yy$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)


imagetitlefig = strcat('results/', fname, '/Ey_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/Ey_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)




%-------------------------

s3=figure()
 

k1 = diff(IMEA3-(IMEA30));
%k1 = diff(IMEA1./(IMEA10+1e-14));

 imagesc(k1)
 
 axis equal;axis tight;
  %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
 ax = gca;
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 
 
 caxis([-1 1]*max(abs(real(k1(:))))) %%%0105 remember to unquote
 % caxis([-0.1 0.1])
 crameri('broc')
 colorbar


Title = strcat('norm g range derivative ($xy$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)


 
imagetitlefig = strcat('results/', fname, '/Exy_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/Exy_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)

 





s3=figure()
 

k1 = diff(IMEA4-(IMEA40));


 imagesc(k1)
 
 axis equal;axis tight;
  %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
 ax = gca;
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 
 
 caxis([-1 1]*max(abs(real(k1(:))))) 
 crameri('broc')
 colorbar


Title = strcat('norm g range derivative ($yx$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)

 
imagetitlefig = strcat('results/', fname, '/Eyx_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/Eyx_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)





%%
%-------------------------REAL G---------------------------------------



s3=figure()
 

k = diff(IMEA1r- IMEA10);
 imagesc(k)
 
 axis equal;axis tight;
 set(gca, 'YAxisLocation', 'left')
 ax = gca;
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 caxis([-1 1]*max(abs(real(k(:))))) 
 crameri('broc')

 colorbar


Title = strcat('norm g (ideal) range derivative ($xx$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title, 'FontSize',13)


 
imagetitlefig = strcat('results/', fname, '/REx_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/REx_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)

 
 
 
s3=figure()


% subplot(211)

 k = diff(IMEA2r-IMEA20);
 %k=imgaussfilt(k,[1.1 1e-12]);
 imagesc(k)
 
 axis equal;axis tight;
 %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
 ax = gca;
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 caxis([-1 1]*max(abs(real(k(:)))))

 crameri('broc')
 colorbar



Title = strcat('norm g (ideal) range derivative ($yy$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)


 
imagetitlefig = strcat('results/', fname, '/REy_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/REy_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)


%---------------

 
s3=figure()


% subplot(211)

 k = diff(IMEA3r-(IMEA30));
 %k=imgaussfilt(k,[1.1 1e-12]);
 imagesc(k)
 
 axis equal;axis tight;
 %set(gca, 'XAxisLocation', 'top')
 set(gca, 'YAxisLocation', 'left')
 ax=gca;
 
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 caxis([-1 1]*max(abs(real(k(:)))))
 crameri('broc')
 colorbar


Title = strcat('norm g (ideal) range derivative ($xy$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)


 
imagetitlefig = strcat('results/', fname, '/RExy_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/RExy_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)






s3=figure()
 

% subplot(211)

k = diff(IMEA4r- IMEA40);
%k = imgaussfilt(k,[1.1 1e-12]);


 imagesc(k)
 
 axis equal;axis tight;
 set(gca, 'YAxisLocation', 'left')
 ax=gca;
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
 
 
 
 
 caxis([-1 1]*max(abs(real(k(:)))))
 crameri('broc')
 colorbar






Title = strcat('norm g (ideal) range derivative ($yx$ component) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',13)


 
imagetitlefig = strcat('results/', fname, '/REyx_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/REyx_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)









%PlotDelta


%% Save image functions
filenamefn = strcat('results/',folder,'/IMGFN_',medium_id,'_',num2str(sep),'_',num2str(tau),'_',num2str(theta));
save(strcat(filenamefn,'_IMEA1.mat'), 'IMEA1');
save(strcat(filenamefn,'_IMEA10.mat'), 'IMEA10');
save(strcat(filenamefn,'_IMEA1r.mat'), 'IMEA1r');
save(strcat(filenamefn,'_IMEA2.mat'), 'IMEA2');
save(strcat(filenamefn,'_IMEA20.mat'), 'IMEA20');
save(strcat(filenamefn,'_IMEA2r.mat'), 'IMEA2r');
save(strcat(filenamefn,'_IMEA3.mat'), 'IMEA3');
save(strcat(filenamefn,'_IMEA30.mat'), 'IMEA30');
save(strcat(filenamefn,'_IMEA3r.mat'), 'IMEA3r');
save(strcat(filenamefn,'_IMEA4.mat'), 'IMEA4');
save(strcat(filenamefn,'_IMEA40.mat'), 'IMEA40');
save(strcat(filenamefn,'_IMEA4r.mat'), 'IMEA4r');






end







