function makeRFfigure(bin_dx_av, bin_resp, figbase)

    %horizontal and vertical components of the dynamic RF
    OF_hor= bin_dx_av(:,:,1);
    OF_ver= bin_dx_av(:,:,2);
    
    %number of rows and columns in the low-res DRF
    [nrow, ncol,xy]=size(bin_dx_av);
    %how many degrees does the half screen span horizontally (used for axis labels)
    half_az_span_deg =70;
    
    
    %create one figure with to plots: Optic flow and strength of RF
    fig=figure('Position',[100,100,1200,500]); 
    %positions of the subplots
    pos1=[0.05,0.15,0.38 0.8];
    pos2=[pos1(1)+pos1(3)+pos1(1)/2,pos1(2),pos1(3), pos1(4)];
    
    %first subplot: OF
    ax1 = subplot('Position',pos1);     
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    lenmax=max(lengths(:));
    OF_hor_plot=OF_hor./lenmax;
    OF_ver_plot=OF_ver./lenmax;
    [x,y] = meshgrid(1:ncol, 1:nrow);
    quiver(x,y,OF_hor_plot,OF_ver_plot, 0,'k','LineWidth',1);     
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);      
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth');
    ylabel('elevation');
    set(ax1,'Position',pos1); 
        
    %second subplot: RF
    ax2 = subplot('Position',pos2);
    bin_resp_plot=bin_resp;
    bin_resp_plot(bin_resp_plot==0)=nan;
    minRF=min(bin_resp(:));
    maxRF=max(bin_resp(:));
    cmap = [1 1 1; parula]; %cells without a value will be white
    imagesc(bin_resp_plot);
    colormap(cmap);
    
    imagesc(bin_resp_plot,[minRF, maxRF]); 
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    colorbar('eastoutside');
    set(ax2, 'Position', pos2)
    xticks([1,ncol/2,ncol]);
    yticks([]);
    yticklabels({});
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth');
    
    ht = suptitle(['RF']);   
    set(ht,'fontsize',12);
    png_name=[figbase,'_RF.png'];    
    saveas(fig,png_name);    
    fig_name=[figbase,'_Rf.fig'];    
    savefig(fig_name);    
end