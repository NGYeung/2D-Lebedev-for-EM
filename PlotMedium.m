function PlotMedium(eps,mu)

figure
subplot(221)
    imagesc(squeeze(eps(:,:,1)))
    axis equal, axis tight
    colorbar
    title('$\varepsilon_{xx}$','Interpreter','latex')
subplot(222)
    imagesc(squeeze(eps(:,:,2)))
    axis equal, axis tight
    colorbar
    title('$\varepsilon_{yy}$','Interpreter','latex')
subplot(223)
    imagesc(squeeze(eps(:,:,3)))
    axis equal, axis tight
    colorbar
    title('$\varepsilon_{xy}$','Interpreter','latex')
subplot(224)
    imagesc(squeeze(mu(:,:)))
    axis equal, axis tight
    colorbar
    title('$\mu$','Interpreter','latex')

end