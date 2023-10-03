%%
 %The delta function
 %plot delta funtions with W in Ex
filtered = 1;


x = 50;
y=140;

filename = strcat('results/',folder,'/DELTAROM_' ,medium_id, '_',num2str(sep), '_', num2str(tau),'_',num2str(x), '_' , num2str(y));



index = sub2ind(sizes.ex1,x,y);





Vex = V(sizes.index_ex1p-sizes.NH,:);
Vey = V(sizes.index_ey2p-sizes.NH,:);
sub = reshape(Vey,[sizes.ey2(1),sizes.ey2(2),size(V,2)]);
sub = sub(:,2:end-1,:);
Vey = reshape(sub,sizes.ex1p,size(V,2));


Vex0 = V0(sizes.index_ex1p-sizes.NH,:);
Vey0 = V0(sizes.index_ey2p-sizes.NH,:);
sub = reshape(Vey0,[sizes.ey2(1),sizes.ey2(2),size(V,2)]);
sub = sub(:,2:end-1,:);
Vey0 = reshape(sub,sizes.ex1p,size(V,2));


Wxx = Vex*Vex0(index,:)';
Wxy = Vex*Vey0(index,:)';
Wyx = Vey*Vex0(index,:)';
Wyy = Vey*Vey0(index,:)';


save(strcat(filename,'xx.mat'), 'Wxx');
save(strcat(filename,'xy.mat'), 'Wxy');
save(strcat(filename,'yx.mat'), 'Wyx');
save(strcat(filename,'yy.mat'), 'Wyy');



clear Vex
clear Vex0
clear Vey
clear Vey0
clear sub

if filtered
    N = 40;
    StartBox=max([x y]-N,1);
    EndBox=min([x y]+N,[sizes.Nx sizes.Ny]);
    

    NxBox=EndBox(1)-StartBox(1);
	NyBox=EndBox(2)-StartBox(2);
	
mask = zeros(sizes.ex1);
mask(StartBox(1):(EndBox(1)-1),(StartBox(2)+1):(EndBox(2)-1)) =1;


index_box = find(mask(:));

 
	
 [msizes, ~,~, meps, ~]=get_medium('Homo',NxBox,NyBox);
 %meps = eps(StartBox(1):(EndBox(1)-1),StartBox(2):(EndBox(2)-1),:);
 [mDEH_sym,mDHE_sym,~]=getSymmetricOperators(msizes,meps);
 AA=mDHE_sym*mDEH_sym;
    AA=-0.5*(AA+AA');
    
    AAalpha_box = AA(msizes.NH + (1:msizes.ex1p),   msizes.NH + (1:msizes.ex1p) );

clear mDEH_sym
clear mDHE_sym
clear AA

[Q,Eig]=eig(full(AAalpha_box),'vector');

fWxx_box=Q * spdiags(sqrt(src.F(sqrt(Eig))), 0, msizes.ex1p, msizes.ex1p) *(Q'*Wxx(index_box));
fWxy_box=Q * spdiags(sqrt(src.F(sqrt(Eig))), 0, msizes.ex1p, msizes.ex1p) *(Q'*Wxy(index_box));
fWyx_box=Q * spdiags(sqrt(src.F(sqrt(Eig))), 0, msizes.ex1p, msizes.ex1p) *(Q'*Wyx(index_box));
fWyy_box=Q * spdiags(sqrt(src.F(sqrt(Eig))), 0, msizes.ex1p, msizes.ex1p) *(Q'*Wyy(index_box));


Wxx_f = zeros(sizes.ex1p,1);Wxy_f = zeros(sizes.ex1p,1);Wyx_f = zeros(sizes.ex1p,1);Wyy_f = zeros(sizes.ex1p,1);

Wxx_f(index_box)=fWxx_box; Wxy_f(index_box)=fWxy_box; Wyx_f(index_box)=fWyx_box; Wyy_f(index_box)=fWyy_box;


%%

f1 = figure()
         Ex1=reshape(Wxx_f,sizes.ex1);
        % Ex1=imgaussfilt(Ex1,[0.1 1e-12]);
         imagesc(real(Ex1));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ex1(:)))))
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
         hold on
         imagesc((eps(:,1:end-1,1)-ones(sizes.ex1)),"AlphaData",0.1)
         crameri('broc')
         axis equal; axis tight;
         title('$\delta_{xx}^{f,ROM}$','FontSize',13)
         colorbar
         ax1 = gca;

f2 = figure()     
         Ey1=reshape(Wxy_f,sizes.ex1);
        % Ey1=imgaussfilt(Ey1,[0.1 1e-12]);
         imagesc(real(Ey1));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ey1(:)))))
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
         hold on
         imagesc(eps(:,1:end-1,3),"AlphaData",0.1)
         crameri('broc')
         axis equal; axis tight;
         title('$\delta_{xy}^{f,ROM}$','FontSize',13)
         colorbar
         ax2 = gca;
         hold off


         
  f3 = figure()          
         Ex2=reshape(Wyx_f,sizes.ex1);
      %   Ex2=imgaussfilt(Ex2,[0.1 1e-12]);
         imagesc(real(Ex2));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ex2(:)))))
         hold on
         imagesc(eps(:,1:end-1,3),"AlphaData",0.1)
         crameri('broc')
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
         axis equal; axis tight;
         title('$\delta_{yx}^{f,ROM}$','FontSize',13)
         colorbar
         ax3 = gca;
         hold off



f4 = figure()          
         Ey2=reshape(Wyy_f,sizes.ex1);
      %   Ey2=imgaussfilt(Ey2,[0.1 1e-12]);
         imagesc(real(Ey2));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ey2(:)))))
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

ylabel(' range in $\lambda$  (x-axis)')
xlabel(' crossrange in $\lambda$  (y-axis)')
         hold on
         imagesc((eps(:,1:end-1,2)-ones(sizes.ex1)),"AlphaData",0.1)
         crameri('broc')
         title('$\delta_{yy}^{f,ROM}$','FontSize',13)
         colorbar
         ax4 = gca;
         hold off


save(strcat(filename,'Filter_xx.mat'), 'Wxx_f');
save(strcat(filename,'Filter_xy.mat'), 'Wxy_f');
save(strcat(filename,'Filter_yx.mat'), 'Wyx_f');
save(strcat(filename,'Filter_yy.mat'), 'Wyy_f');



saveas(f1,strcat(filename,'xx_filtered.fig'))
saveas(f2,strcat(filename,'xy_filtered.fig'))
saveas(f3,strcat(filename,'yx_filtered.fig'))
saveas(f4,strcat(filename,'yy_filtered.fig'))
exportgraphics(ax1,strcat(filename,'xx_filtered.png'),'Resolution',300)
exportgraphics(ax2,strcat(filename,'xy_filtered.png'),'Resolution',300)
exportgraphics(ax3,strcat(filename,'yx_filtered.png'),'Resolution',300)
exportgraphics(ax4,strcat(filename,'yy_filtered.png'),'Resolution',300)

return
end








%%
% unfiltered

f1 = figure()
         Ex1=reshape(Wxx,sizes.ex1);
         imagesc(real(Ex1));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ex1(:)))))
         hold on
         imagesc((eps(:,1:end-1,1)-ones(sizes.ex1)),"AlphaData",0.1)
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})
         crameri('broc')
         axis equal; axis tight;
         title('$\delta_{xx}^{ROM}$','FontSize',13)
         colorbar
         ax1=gca;
         hold off

f2 = figure()     
         Ey1=reshape(Wxy,sizes.ex1);
         imagesc(real(Ey1));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ey1(:)))))
         hold on
         imagesc(eps(:,1:end-1,3),"AlphaData",0.1)
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})
         crameri('broc')
         axis equal; axis tight;
         title('$\delta_{xy}^{ROM}$','FontSize',13)
         colorbar
         ax2=gca;
         hold off
         
  f3 = figure()          

         Ex2=reshape(Wyx,sizes.ex1);
         imagesc(real(Ex2));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ex2(:)))))
         hold on
         imagesc(eps(:,1:end-1,3),"AlphaData",0.1)
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})
         crameri('broc')
         axis equal; axis tight;
         title('$\delta_{yx}^{ROM}$','FontSize',13)
         colorbar
         ax3=gca;
         hold off

f4 = figure()          

         Ey2=reshape(Wyy,sizes.ex1);
         imagesc(real(Ey2));
         crameri('broc')
         axis equal; axis tight;
         caxis([-1 1]*max(abs(real(Ey2(:)))))
         hold on
         imagesc((eps(:,1:end-1,2)-ones(sizes.ex1)),"AlphaData",0.1)
         xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})
         crameri('broc')
         title('$\delta_{yy}^{ROM}$','FontSize',13)
         colorbar
         ax4=gca;
         hold off


saveas(f1,strcat(filename,'xx.fig'))
saveas(f2,strcat(filename,'xy.fig'))
saveas(f3,strcat(filename,'yx.fig'))
saveas(f4,strcat(filename,'yy.fig'))
exportgraphics(ax1,strcat(filename,'xx.png'),'Resolution',300)
exportgraphics(ax2,strcat(filename,'xy.png'),'Resolution',300)
exportgraphics(ax3,strcat(filename,'yx.png'),'Resolution',300)
exportgraphics(ax4,strcat(filename,'yy.png'),'Resolution',300)
% 


return

