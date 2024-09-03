function make_figure_RF_one_cell(cell_data,figname)    
    xvals = cell_data.maxVlines;
    yvals = cell_data.maxHlines;
    num_hor_bins = cell_data.maxVlines;
    num_ver_bins = cell_data.maxHlines;
    
    RFx_sub=imresize(cell_data.RF_x,[num_ver_bins,num_hor_bins],"nearest");
    RFy_sub=imresize(cell_data.RF_y,[num_ver_bins,num_hor_bins],"nearest");
    

    RF_amp=sqrt(RFx_sub.^2.+RFy_sub.^2);
    lenmax=max(RF_amp(:));
    OF_hor_plot=RFx_sub./lenmax;
    OF_ver_plot=RFy_sub./lenmax;


   f = figure; t=tiledlayout(2,2,"TileSpacing","compact");
    %how many degrees does the half screen span horizontally (used for axis labels)
    half_az_span_deg =70; 
    nexttile;
    [x,y] = meshgrid(1:num_hor_bins,1:num_ver_bins);
    %contour(RF_patch,[0.5 0.5],'-b'); hold on;    
    quiver(x,y,OF_hor_plot,OF_ver_plot, 0,'k','LineWidth',1);     
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 num_hor_bins+1]);
    ylim([0 num_ver_bins+1]);      
    xticks([1,num_hor_bins/2,num_hor_bins]);
    yticks([num_ver_bins/2]);
    yticklabels(0);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth');
    ylabel('elevation');
    
    nexttile;  
    maxRF=max(abs(RF_amp(:)));
    minRF=-maxRF;
    if any(isnan(RF_amp(:)))
        cmap = [1 1 1; parula]; %cells without a value will be white
    else
        cmap = parula;
    end    
    colormap(cmap);
    if any(isnan([minRF, maxRF]))
        imagesc(RF_amp);  
    else
        imagesc(RF_amp,[minRF, maxRF]); 
    end
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    axis equal;
    xlim([0 num_hor_bins+1]);
    ylim([0 num_ver_bins+1]);
    colorbar('eastoutside');    
    xticks([1,num_hor_bins/2,num_hor_bins]);
    yticks([]);
    yticklabels({});
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    xlabel('azimuth'); 

    nexttile;
    maxRF=max(abs(RFx_sub(:)));
    minRF=-maxRF;
    if any(isnan(RFx_sub(:)))
        cmap = [1 1 1; parula]; %cells without a value will be white
    else
        cmap = parula;
    end    
    colormap(cmap);
    if any(isnan([minRF, maxRF]))
        imagesc(RFx_sub);  
    else
        imagesc(RFx_sub,[minRF, maxRF]); 
    end
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    axis equal;
    xlim([0 num_hor_bins+1]);
    ylim([0 num_ver_bins+1]);
    colorbar('eastoutside');    
    xticks([1,num_hor_bins/2,num_hor_bins]);
    yticks([]);
    yticklabels({});
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);

    nexttile;
    maxRF=max(abs(RFy_sub(:)));
    minRF=-maxRF;
    if any(isnan(RFy_sub(:)))
        cmap = [1 1 1; parula]; %cells without a value will be white
    else
        cmap = parula;
    end    
    colormap(cmap);
    if any(isnan([minRF, maxRF]))
        imagesc(RFy_sub);  
    else
        imagesc(RFy_sub,[minRF, maxRF]); 
    end
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    axis equal;
    xlim([0 num_hor_bins+1]);
    ylim([0 num_ver_bins+1]);
    colorbar('eastoutside');    
    xticks([1,num_hor_bins/2,num_hor_bins]);
    yticks([]);
    yticklabels({});
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);    

    titlestr=[cell_data.strain,' ', cell_data.celltype,' ',cell_data.datestr, ' ', cell_data.cellid];
    title(t,titlestr, 'Interpreter', 'none');
    if ~isempty(figname)       
        savefig(figname);
    end
    close(f);
end