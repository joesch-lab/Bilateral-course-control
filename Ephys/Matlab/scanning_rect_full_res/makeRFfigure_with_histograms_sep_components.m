function makeRFfigure_with_histograms_sep_components(bin_dx_av, bin_resp, bin_hor, bin_ver, figbase)
    %horizontal and vertical components of the dynamic RF
    OF_hor= bin_dx_av(:,:,1);
    OF_ver= bin_dx_av(:,:,2);
    
    %number of rows and columns in the low-res DRF
    [nrow, ncol,xy]=size(bin_dx_av);
    %how many degrees does the half screen span horizontally (used for axis labels)
    half_az_span_deg =70;
    
    
    %create one figure with to plots: Optic flow and strength of RF
    fig=figure('Position',[100,100,800,600]); 
    
    %positions of the subplots
    %margin btw subplots
    marg=0.05;
    %width and height of the OF subplot
    plotthick=0.7-2*marg;
    plotthick_min=plotthick;
    %height of the bar subplot
    barthick=0.25-marg;
    
    %positions of subplots    
    pos1=[marg,2*marg+barthick,barthick, plotthick_min];
    pos2=[2*marg+barthick,marg,plotthick,barthick];
    pos3=[2*marg+barthick,2*marg+barthick,plotthick,plotthick_min];
    
    %horizontal bars
    ax1 = subplot('Position',pos1);     
    nvb=length(bin_ver);
    barh(1:nvb,bin_ver');
    xlim([min(bin_ver(:)),max(bin_ver(:))]);
    xticks([min(bin_ver(:)),max(bin_ver(:))]);
    yticks([]);    
    ylabel('elevation');
       
    
    %subplot: OF
    ax3 = subplot('Position',pos3);     
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    lenmax=max(lengths(:));
    OF_hor_plot=OF_hor./lenmax;
    OF_ver_plot=OF_ver./lenmax;
    [x,y] = meshgrid(1:ncol, 1:nrow);
    quiver(x,y,OF_hor_plot,OF_ver_plot, 0,'k','LineWidth',1);     
    set(gca,'XDir','reverse');     
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);      
    axis equal;
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    linkaxes([ax1,ax3],'y');
    
    %vertical bars
    ax2 = subplot('Position',pos2);     
    nhb=length(bin_hor);
    bar(1:nhb,bin_hor');
    ax2.YAxisLocation = 'right';
    xticks([]);    
    ylim([min(bin_hor(:)),max(bin_hor(:))]);
    yticks([min(bin_hor(:)),max(bin_hor(:))]);
    set(gca,'XDir','reverse'); 
    xlabel('azimuth');
    linkaxes([ax3,ax2],'x');
    
    png_name=[figbase,'_RF_with_distr.png'];    
    saveas(fig,png_name);    
    fig_name=[figbase,'_Rf.fig'];    
    savefig(fig_name);
end