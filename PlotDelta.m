%%
 %The delta function
 %plot delta funtions with W in Ex
filtered = 0;


x = 95;
y = 130;
filename = strcat('results/',folder,'/DELTAROM_' ,medium_id, '_',num2str(sep), '_', num2str(tau),'_',num2str(x), '_' , num2str(y));



index = sub2ind(sizes.ex1,x,y);



if filtered

A_f = zeros(sizes.ex1p+sizes.ey1p);
A_f(1:sizes.ex1p,1:sizes.ex1p) = Ac(sizes.index_ex1p-sizes.NH,sizes.index_ex1p-sizes.NH);
A_f(sizes.ex1p+(1:sizes.ey2p),1:sizes.ex1p) = Ac(sizes.index_ey2p-sizes.NH,sizes.index_ex1p-sizes.NH);
A_f(1:sizes.ex1p,sizes.ex1p+(1:sizes.ey2p)) = Ac(sizes.index_ex1p-sizes.NH,sizes.index_ey2p-sizes.NH);
A_f(sizes.ex1p+(1:sizes.ey2p),sizes.ex1p+(1:sizes.ey2p)) = Ac(sizes.index_ey2p-sizes.NH,sizes.index_ey2p-sizes.NH);

A_f = sparse(A_f);
disp('Compute the decomposition for delta function.')


tic;[Q,lambda0ftr]=eig(full(A_f),'vector');toc
bSmall=Q * spdiags(sqrt(src.F(sqrt(lambda0ftr))), 0, 2*sizes.ex1p, 2*sizes.ex1p) *Q';
clear A_f



V_F = bSmall*V;

Vfex = V_F(sizes.index_ex1p-sizes.NH,:);
Vfey = V_F(sizes.index_ey2p-sizes.NH,:);
sub = reshape(Vfey,[sizes.ey2(1),sizes.ey2(2),size(V,2)]);
sub = sub(:,2:end-1,:);
Vfey = reshape(sub,sizes.ex1p,size(V,2));


clear V_F
clear bSmall


end




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



if filtered
Wxx_f = Vfex*Vex0(index,:)';
Wxy_f = Vfex*Vey0(index,:)';
Wyx_f = Vfey*Vex0(index,:)';
Wyy_f = Vfey*Vey0(index,:)';
end


clear Vex
clear Vex0
clear Vey
clear Vey0
clear sub

if filtered
clear Vfey 
clear Vfex
end




%%
%filtered
if filtered

f1 = figure()
         Ex1=reshape(Wxx_f,sizes.ex1);
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



saveas(f1,strcat(filename,'xx_filtered.fig'))
saveas(f2,strcat(filename,'xy_filtered.fig'))
saveas(f3,strcat(filename,'yx_filtered.fig'))
saveas(f4,strcat(filename,'yy_filtered.fig'))
exportgraphics(ax1,strcat(filename,'xx_filtered.fig'),'Resolution',300)
exportgraphics(ax2,strcat(filename,'xy_filtered.fig'),'Resolution',300)
exportgraphics(ax3,strcat(filename,'yx_filtered.fig'),'Resolution',300)
exportgraphics(ax4,strcat(filename,'yy_filtered.fig'),'Resolution',300)

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

