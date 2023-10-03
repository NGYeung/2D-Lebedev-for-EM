%plot sum
s3=figure()




 k = diff(IMEA1-(IMEA10)+IMEA4-(IMEA40));

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


Title = strcat('norm g range derivative ( from $x$ polarized source) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',12)


imagetitlefig = strcat('results/', fname, '/X_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/X_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)




%-------------------------

s3=figure()
 

k1 = diff(IMEA3-(IMEA30)+IMEA2-(IMEA20));
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


Title = strcat('norm g range derivative (from $y$ polarized sources) ', '$\tau$ = ', num2str(r),'$\pi/w_c$');
title(Title,'FontSize',12)


 
imagetitlefig = strcat('results/', fname, '/RY_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.fig');
imagetitlepng = strcat('results/', fname, '/RY_image_',  medium_id, '_diff_', num2str(tau), '_', num2str(sep),'_', num2str(theta), '.png');
saveas(s3,imagetitlefig)
exportgraphics(ax,imagetitlepng,'Resolution',300)
