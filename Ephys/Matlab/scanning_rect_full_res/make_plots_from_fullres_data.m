%use this script to load full res data and make plots with the given number of bins

%this file should contain the vaules for
% 'V_frame_resp','bin_edges','Data', 'stim_arr', 'bin_dx_av','bin_resp',
% 'RFx_hres', 'RFy_hres', 'RFx_pos_hres', 'RFx_neg_hres', 'RFy_pos_hres', 'RFy_neg_hres'

% datamatfile='C:\DATA\Data_for_Olga\scanning_rect\stimuli_scanning_rect_rescale_with_pause_2020_09_30_15_02_49.mat';
% load(datamatfile);
% 
% %specify how many horizontal bins in the new plot
% num_hor_bins=13;
% [bin_dx_av, bin_resp] = deform_plot_uniform_sampling(RFx_hres, RFy_hres, num_hor_bins, 'c', 0);
%   
% [filepath,name,ext] = fileparts(datamatfile);    
% figbase = fullfile(filepath,name);
makeRFfigure(bin_dx_av, bin_resp, figbase);
makeRF_indiv_compfigure(RFx_pos_hres, RFx_neg_hres, RFy_pos_hres, RFy_neg_hres, num_hor_bins, figbase);

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

    
function makeRF_indiv_compfigure(RFx_pos_hres, RFx_neg_hres, RFy_pos_hres, RFy_neg_hres, num_hor_bins, figbase)
    
    emptyRF = zeros(size(RFx_pos_hres));
    [pos_dx_av, pos_dx_resp] = deform_plot_uniform_sampling(RFx_pos_hres, emptyRF, num_hor_bins, 'x', 0);
    [neg_dx_av, neg_dx_resp] = deform_plot_uniform_sampling(RFx_neg_hres, emptyRF, num_hor_bins, 'x', 0);
    [pos_dy_av, pos_dy_resp] = deform_plot_uniform_sampling(emptyRF,   RFy_pos_hres, num_hor_bins, 'y', 0);
    [neg_dy_av, neg_dy_resp] = deform_plot_uniform_sampling(emptyRF,   RFy_neg_hres, num_hor_bins, 'y', 0);
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
    
    subplot_tight(2,4,1,[hm1,vm1]);  
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
    
    subplot_tight(2,4,2,[hm1,vm1]);   
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
    
    subplot_tight(2,4,3,[hm1,vm1]);   
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
    
    subplot_tight(2,4,4,[hm1,vm1]);   
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
    subplot_tight(2,4,5,[hm2,vm1]); 
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
    
    subplot_tight(2,4,6,[hm2,vm1]);    
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
    
    subplot_tight(2,4,7,[hm2,vm1]); 
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
    
    subplot_tight(2,4,8,[hm2,vm1]);       
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

