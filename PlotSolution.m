%This script plots the internal solution


index = [60,210]; %change the index.
%index = [30,30];

maskx = zeros(sizes.ex1(1),sizes.ex1(2));
maskx(index(1),index(2))=1;
indexx = find(maskx(:));

masky = zeros(sizes.ey2(1),sizes.ey2(2));
masky(index(1),index(2))=1;
indexy = sizes.ex1p + find(masky(:));

%% x-polarized field
Solx = V0(indexx,:)*R;
%Solx = U_E(indexx,:)-U_E0(indexx,:);
Solx = reshape(Solx,src.nsrc,size(D,3)/2);

%% y-polarized solution
Soly = V0(indexy,:)*R
%Soly = U_E(indexy,:)-U_E0(indexy,:);
Soly = reshape(Soly,src.nsrc,size(D,3)/2);

filenamesol = strcat('results/',folder,'/IMGFN_',medium_id,'_',num2str(sep),'_',num2str(tau),'_',num2str(x),'_',num2str(y));
save(strcat(filenamesol,'_solx.mat'), 'Solx');
save(strcat(filenamesol,'_soly.mat'), 'Soly');


%%

s6 =figure(31)


subplot(121)
imagesc(Solx)

hold on
plot([0 72],[80 80])

xlabel('time steps')
ylabel('sensor index. 1-80: x-polarized; 81-160: y-polarized ')

caxis([-1 1]*max(abs(real(Solx(:)))))
crameri('broc')
title('Internal wave (x-polarized). ')
colorbar



subplot(122)
imagesc(Soly)

hold on
plot([0 72],[80 80])


xlabel('time steps')
ylabel('sensor index. 1-80: x-polarized; 81-160: y-polarized ')
 
caxis([-1 1]*max(abs(real(Soly(:)))))
crameri('broc')
title('Internal wave  (y-polarized)')
colorbar

ax = gcf;

Title = strcat('Internal Wave $g(jt,x_r;y)$ at point (',num2str(index(1)),', ',num2str(index(2)),')' );
sgtitle(Title,'FontSize',14)

imagetitlefigs = strcat('results/', fname, '/Internal_', medium_id, '_nbk_' , num2str(tau), '.fig');
imagetitlepngs = strcat('results/', fname, '/Internal_', medium_id, '_nbk_' , num2str(tau), '.png');
saveas(s6,imagetitlefigs) 
exportgraphics(ax,imagetitlepngs)





