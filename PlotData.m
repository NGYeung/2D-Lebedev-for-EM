%D_j at time j is calculated as b^Tu_j where u_j is Ntot x source (2m x Ntot)*(Ntot x 1)

%source_i=2;
source_i = 40;
% source/receiver # and polarization. 1 = x, 2 = y; %central sensor
datas = zeros(src.nsrc,size(D,3));
datasog = zeros(src.nsrc,size(D,3));


filenamefig = strcat('results/', folder, '/' ,medium_id, '_data0_',num2str(source_i) , '_', num2str(Nx*Ny), '_', num2str(tau), '.fig');
filenamepng = strcat('results/', folder, '/' ,medium_id, '_data0_',num2str(source_i) , '_', num2str(Nx*Ny), '_', num2str(tau), '.png');

filenamefigm = strcat('results/', folder, '/' ,medium_id, '_data_',num2str(source_i) , '_', num2str(Nx*Ny), '_', num2str(tau), '.fig');
filenamepngm = strcat('results/', folder, '/' ,medium_id, '_data_',num2str(source_i) , '_', num2str(Nx*Ny), ' ', num2str(tau), '.png');



for time = 1:size(D,3)                             
datas(:,time) = D(:,source_i,time)-D0(:,source_i,time);
datasog(:,time) = D(:,source_i,time);
end

s5 = figure()
imagesc(real(datas))

hold on
plot([0 144],[80 80])

axis equal;axis tight;
xlabel('Time steps j')
ylabel('sensor index')
caxis([-1 1].*max(abs(real(datas(:)))))
crameri('broc')
title('$D_j-D_{0,j}$ from central x-polarized sensor','FontSize',13)
colorbar
ax = gca;
saveas(s5,filenamefigm)
exportgraphics(ax,filenamepngm,'Resolution',300)





s5 = figure()
imagesc(real(datasog))
hold on
plot([0 144],[80 80])

axis equal;axis tight;
xlabel('Time steps j')
ylabel('sensor index')
caxis([-1 1].*max(abs(real(datasog(:)))))
crameri('broc')
title('$D_j$ from central x-polarized sensor','FontSize',13)
colorbar
ax = gca;
saveas(s5,filenamefig)
exportgraphics(ax,filenamepng,'Resolution',300)