function makeRF_indiv_compfigure(RFx_pos_hres, RFx_neg_hres, RFy_pos_hres, RFy_neg_hres, num_hor_bins, figbase)
    
    emptyRF = zeros(size(RFx_pos_hres));
    [pos_dx_av, pos_dx_resp, pos_dx_resp_deformed] = deform_plot_uniform_sampling(RFx_pos_hres, emptyRF, num_hor_bins, 'x', 0);
    [neg_dx_av, neg_dx_resp, neg_dx_resp_deformed] = deform_plot_uniform_sampling(RFx_neg_hres, emptyRF, num_hor_bins, 'x', 0);
    [pos_dy_av, pos_dy_resp, pos_dy_resp_deformed] = deform_plot_uniform_sampling(emptyRF,   RFy_pos_hres, num_hor_bins, 'y', 0);
    [neg_dy_av, neg_dy_resp, neg_dy_resp_deformed] = deform_plot_uniform_sampling(emptyRF,   RFy_neg_hres, num_hor_bins, 'y', 0);
    %set to nan zero values of the RF    
    pos_dx_resp(pos_dx_resp==0)=nan; neg_dx_resp(neg_dx_resp==0)=nan; 
    pos_dy_resp(pos_dy_resp==0)=nan; neg_dy_resp(neg_dy_resp==0)=nan;
    
    [nrow, ncol,xy]=size(pos_dx_av);
    half_az_span_deg =70;    
    
    %find the max length of the vector to normalize all the directions to
    %the unit length
    maxlen=0;
    OF_hor= pos_dx_av(:,:,1); OF_ver= pos_dx_av(:,:,2);
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    maxlen=max(maxlen, max(lengths(:)));
    OF_hor= neg_dx_av(:,:,1); OF_ver= neg_dx_av(:,:,2);
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    maxlen=max(maxlen, max(lengths(:)));
    OF_hor= pos_dy_av(:,:,1); OF_ver= pos_dy_av(:,:,2);
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    maxlen=max(maxlen, max(lengths(:)));
    OF_hor= neg_dy_av(:,:,1); OF_ver= neg_dy_av(:,:,2);
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    maxlen=max(maxlen, max(lengths(:)));
    
    %rescale responses according to the deformation: responses on the sides
    %where xy-projection of the stimulus is small, get smaller weights
    pos_dx_resp=pos_dx_resp.*(pos_dx_resp_deformed/max(pos_dx_resp_deformed(:)));
    neg_dx_resp=neg_dx_resp.*(neg_dx_resp_deformed/max(neg_dx_resp_deformed(:)));
    pos_dy_resp=pos_dy_resp.*(pos_dy_resp_deformed/max(pos_dy_resp_deformed(:)));
    neg_dy_resp=neg_dy_resp.*(neg_dy_resp_deformed/max(neg_dy_resp_deformed(:)));
    
    %min and max responses
    min_resp = min([min(pos_dx_resp(:)), min(neg_dx_resp(:)), min(pos_dy_resp(:)), min(neg_dy_resp(:))]);
    max_resp = max([max(pos_dx_resp(:)), max(neg_dx_resp(:)), max(pos_dy_resp(:)), max(neg_dy_resp(:))]);

    %empty 
    hones=zeros(size(pos_dx_resp));
    [x,y] = meshgrid(1:ncol, 1:nrow);
    %colomap
    cmap = [1 1 1; parula]; %cells without a value will be white    
    
    fig=figure('Position',[100,100,1200,600]); 
    vm1=0.04;
    hm1=0.001;
    marg=0.05;
    topmarg=0.1;
    plotw=(1-marg)/4-marg;
    ploth=(1-topmarg)/2-marg;
    
    pos1=[marg,marg,plotw, ploth];
    pos2=[2*marg+plotw,marg,plotw, ploth];
    pos3=[3*marg+2*plotw,marg,plotw, ploth];
    pos4=[4*marg+3*plotw,marg,plotw, ploth];
    
    pos5=[marg,2*marg+ploth,plotw, ploth];
    pos6=[2*marg+plotw,2*marg+ploth,plotw, ploth];
    pos7=[3*marg+2*plotw,2*marg+ploth,plotw, ploth];
    pos8=[4*marg+3*plotw,2*marg+ploth,plotw, ploth];
    
    
    
%     subplot_tight(2,4,1,[hm1,vm1]);  
    ax5 = subplot('Position',pos5);    
    quiver(x,y,pos_dx_av(:,:,1)./maxlen,hones, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    ylabel('elevation');    
    title('Leftward')
    
%     subplot_tight(2,4,2,[hm1,vm1]);  
    ax6 = subplot('Position',pos6);    
    quiver(x,y,-neg_dx_av(:,:,1)./maxlen,hones, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    title('Rightward');
    
%     subplot_tight(2,4,3,[hm1,vm1]);
    ax7 = subplot('Position',pos7);
    quiver(x,y, hones, pos_dy_av(:,:,2)./maxlen, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    title('Upward');
    
%     subplot_tight(2,4,4,[hm1,vm1]);   
    ax8 = subplot('Position',pos8);
    quiver(x,y, hones, -neg_dy_av(:,:,2)./maxlen, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
%     xlabel('azimuth'); 
    title('Downward');    
    
    %RF strength
    hm2=0.001;
%     subplot_tight(2,4,5,[hm2,vm1]); 
    ax1 = subplot('Position',pos1);
    imagesc(pos_dx_resp,[min_resp, max_resp]); 
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth'); 
    ylabel('elevation');  
    axis equal;
    axis tight;   
    
%     subplot_tight(2,4,6,[hm2,vm1]); 
    ax2 = subplot('Position',pos2);
    imagesc(neg_dx_resp,[min_resp, max_resp]); 
    set(gca,'YDir','normal');  
    set(gca,'Xdir','reverse');
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth');     
    axis equal;
    axis tight;
    
%     subplot_tight(2,4,7,[hm2,vm1]); 
    ax3 = subplot('Position',pos3);
    imagesc(pos_dy_resp,[min_resp, max_resp]); %need to flip left right for the fly's space
    set(gca,'YDir','normal');    
    set(gca,'Xdir','reverse');
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth');  
    axis equal;
    axis tight;
    
%     subplot_tight(2,4,8,[hm2,vm1]);   
    ax4 = subplot('Position',pos4);
    imagesc(neg_dy_resp,[min_resp, max_resp]); %need to flip left right for the fly's space
    set(gca,'Xdir','reverse');
    set(gca,'YDir','normal');
    pos1=get(gca,'Position');
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    xticks([1,ncol/2,ncol]);
    yticks([nrow/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth');  
    axis equal;
    axis tight;
    colormap(cmap);
   
    cbar_handle =  colorbar('westoutside');   
    set(cbar_handle, 'YAxisLocation','right'); 
    pos21=pos1(1)+pos1(3)+0.01;
    pos22=0.125;%pos1(2)*10;
    pos23=0.01;
    pos24=pos1(4)/2;
    set(cbar_handle,'position',[pos21, pos22,pos23,pos24]);
    suptitle(['Individual components of the dynamic RF']); 
    RF_name=[figbase,'_ind_components_RFOF.png'];
    saveas(fig,RF_name); 
    fig_name=[figbase,'_ind_components_RFOF.fig'];    
    savefig(fig_name);       
end