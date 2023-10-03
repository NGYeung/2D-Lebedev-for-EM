s2 = figure(21) 
%missing source location here. fix tmr.

        subplot(222)
        hold on
        
        Ey=eps(:,:,2);
        
        imagesc((real(Ey)));
        colormap(flipud(pink))
        plot(src.src_loc(1,:), src.src_loc(2,:), 'x')
        %plot(80,90,'+')

        set(gca,'YDir','reverse')
        set(gca,'XDir','normal')
        axis equal;axis tight
        caxis([0 1]*max(abs(real(Ey(:)))))
        
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})
 
ylabel(' range in $\lambda$ ')
xlabel(' crossrange in $\lambda$ ')
 
        rectangle('Position',[20 20 sizes.Ny-40 sizes.Nx-60]);
        title('epsilon yy')
        colorbar
        hold off

        subplot(221)
        hold on
        Ex=eps(:,:,1);
        imagesc((real(Ex)));
        colormap(flipud(pink))
        plot(src.src_loc(1,:), src.src_loc(2,:), 'x')
        %plot(80,90,'+')

        set(gca,'YDir','reverse')
        set(gca,'XDir','normal')
        
        axis equal;axis tight
        caxis([0 1]*max(abs(real(Ex(:)))))
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

 
ylabel(' range in $\lambda$ ')
xlabel(' crossrange in $\lambda$ ')
 
        rectangle('Position',[20 20 sizes.Ny-40 sizes.Nx-60]);
        title('epsilon xx')
        colorbar
        hold off
        
        
        subplot(223)
        hold on
        Exy=eps(:,:,3);
        imagesc((real(Exy)));
        colormap(flipud(pink))
        plot(src.src_loc(1,:), src.src_loc(2,:), 'x')
        %plot(80,90,'+')

        set(gca,'YDir','reverse')
        set(gca,'XDir','normal')
        axis equal;axis tight
        caxis([0 1]*max(abs(real(Ex(:)))))
xticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
xticklabels({'0','2', '4', '6' , '8', '10', '12'})

yticks([0  2*lamb_c  4*lamb_c  6*lamb_c 8*lamb_c 10*lamb_c  12*lamb_c ])
yticklabels({'0','2', '4', '6' , '8', '10', '12'})

 
ylabel(' range in $\lambda$ ')
xlabel(' crossrange in $\lambda$ ')
 
        rectangle('Position',[20 20 sizes.Ny-40 sizes.Nx-60]);
        title('epsilon xy')
        colorbar
        hold off
        
        
      
        
        sgtitle('Configuration')
        
saveas(s2,strcat('results/', fname, '/config_image_',medium_id, '.fig'))
ax = gcf;
exportgraphics(ax, strcat('results/', fname, '/config_image_',medium_id, '.png'),'Resolution',300)


