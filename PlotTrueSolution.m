%This script plots the internal solution


x =60;
y = 210;

index = [x,y]; %change the index.
%index = [30,30];

maskx = zeros(sizes.ex1(1),sizes.ex1(2));
maskx(index(1),index(2))=1;
indexx = find(maskx(:));

masky = zeros(sizes.ey2(1),sizes.ey2(2));
masky(index(1),index(2))=1;
indexy = sizes.ex1p + find(masky(:));

%% x-polarized field
%Solx = V0(indexx,:)*R;
Solx = U_E(indexx,:);
Solx0 = U_E0(indexx,:);
Solx = reshape(Solx,src.nsrc,size(D,3)/2);
Solx0 = reshape(Solx0,src.nsrc,size(D,3)/2);

%% y-polarized solution
%Soly = V0(indexy,:)*R
Soly = U_E(indexy,:);
Soly0 = U_E0(indexy,:);
Soly = reshape(Soly,src.nsrc,size(D,3)/2);
Soly0 = reshape(Soly0,src.nsrc,size(D,3)/2);

filenamesol = strcat('results/',folder,'/Internal_',medium_id,'_',num2str(sep),'_',num2str(tau),'_',num2str(x),'_',num2str(y));
save(strcat(filenamesol,'_solx.mat'), 'Solx');
save(strcat(filenamesol,'_soly.mat'), 'Soly');
save(strcat(filenamesol,'_solx0.mat'), 'Solx0');
save(strcat(filenamesol,'_soly0.mat'), 'Soly0');

%%

s6 =figure(32)


subplot(121)
imagesc(Solx)

hold on
plot([0 72],[80 80])

xlabel('time steps')
ylabel('sensor index. 1-80: x-polarized; 81-160: y-polarized ')

caxis([-1 1]*max(abs(real(Solx(:)))))
crameri('broc')
title('True Internal wave (x-polarized). ')
colorbar



subplot(122)
imagesc(Soly)

hold on
plot([0 72],[80 80])


xlabel('time steps')
ylabel('sensor index. 1-80: x-polarized; 81-160: y-polarized ')
 
caxis([-1 1]*max(abs(real(Soly(:)))))
crameri('broc')
title('True Internal wave  (y-polarized)')
colorbar

ax = gcf;

Title = strcat('Internal Wave $g(jt,x_r;y)$ at point (',num2str(index(1)),', ',num2str(index(2)),')' );
sgtitle(Title,'FontSize',14)

imagetitlefigs = strcat('results/', fname, '/Internal_true_', medium_id, num2str(tau), '.fig');
imagetitlepngs = strcat('results/', fname, '/Internal_true_', medium_id, num2str(tau), '.png');
saveas(s6,imagetitlefigs) 
exportgraphics(ax,imagetitlepngs)




%%
%--------------------------------------------------------------------------------------

s6 =figure(31)


subplot(121)
imagesc(Solx0)

hold on
plot([0 72],[80 80])

xlabel('time steps')
ylabel('sensor index. 1-80: x-polarized; 81-160: y-polarized ')

caxis([-1 1]*max(abs(real(Solx0(:)))))
crameri('broc')
title('Internal wave (background, x-polarized). ')
colorbar



subplot(122)
imagesc(Soly0)

hold on
plot([0 72],[80 80])


xlabel('time steps')
ylabel('sensor index. 1-80: x-polarized; 81-160: y-polarized ')
 
caxis([-1 1]*max(abs(real(Soly0(:)))))
crameri('broc')
title('Internal wave  (background, y-polarized)')
colorbar

ax = gcf;

Title = strcat('Internal Wave $g(jt,x_r;y)$ at point (',num2str(index(1)),', ',num2str(index(2)),')' );
sgtitle(Title,'FontSize',14)

imagetitlefigs = strcat('results/', fname, '/Internal_bkg_', medium_id, num2str(tau), '.fig');
imagetitlepngs = strcat('results/', fname, '/Internal_bkg_', medium_id, num2str(tau), '.png');
saveas(s6,imagetitlefigs) 
exportgraphics(ax,imagetitlepngs)



return
