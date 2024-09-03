function make_figure_RF_one_rep(rep_data,resfolder, extra_remark)
    xvals = rep_data.stim_info.scan_v_x;
    yvals = rep_data.stim_info.scan_h_y;
    num_hor_bins = length(xvals);
    num_ver_bins = length(yvals);
    
    RFx_sub=zeros(num_ver_bins,num_hor_bins);
    RFy_sub=zeros(num_ver_bins,num_hor_bins);
    for yi=1:num_ver_bins 
        RFx_sub(yi,:)=interp1(rep_data.stim_info.scan_h_x,rep_data.RF.RFx(yi,:),xvals,"nearest");
    end
    for xi=1:num_hor_bins 
        RFy_sub(:,xi)=interp1(rep_data.stim_info.scan_v_y,rep_data.RF.RFy(:,xi),yvals,"nearest");
    end
    


%     [nrow,ncol]=size(rep_data.RF.RF);
%     num_ver_bins=ceil(nrow*num_hor_bins/ncol);
%     hor_edges=round(linspace(1,ncol,num_hor_bins+1));
%     ver_edges=round(linspace(1,nrow,num_ver_bins+1));
% 
%     RF_lowres=zeros(num_ver_bins,num_ver_bins,3);
%     %compute mean x,y responses and patch values
%     for ri=1:num_ver_bins
%         for ci=1:num_hor_bins
%             subRF = rep_data.RF.RFx(ver_edges(ri):ver_edges(ri+1),hor_edges(ci):hor_edges(ci+1));
%             RF_lowres(ri,ci,1)=mean(subRF(:));
%             subRF = rep_data.RF.RFy(ver_edges(ri):ver_edges(ri+1),hor_edges(ci):hor_edges(ci+1));
%             RF_lowres(ri,ci,2)=mean(subRF(:));            
%         end
%     end
% 
%     lengths=sqrt(RF_lowres(:,:,1).^2.+RF_lowres(:,:,2).^2);
%     lenmax=max(lengths(:));
%     OF_hor_plot=RF_lowres(:,:,1)./lenmax;
%     OF_ver_plot=RF_lowres(:,:,2)./lenmax;
% 
% 

    RF_amp=sqrt(RFx_sub.^2.+RFy_sub.^2);
    lenmax=max(RF_amp(:));
    OF_hor_plot=RFx_sub./lenmax;
    OF_ver_plot=RFy_sub./lenmax;


   f = figure; t=tiledlayout(3,3,"TileSpacing","compact");
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
   

    nexttile; axis('off');

    nexttile(4,[1,2]);
    xvals=(1:length(rep_data.trace.htrace))/10000;
    maxval=max(rep_data.trace.htrace);    
    lefton(rep_data.trace.h_left_t)=maxval;
    lefton=zeros(1,length(rep_data.trace.htrace));
    lefton(rep_data.trace.h_left_t)=maxval;
    righton=zeros(1,length(rep_data.trace.htrace));
    righton(rep_data.trace.h_right_t)=maxval;
    redsignal=rep_data.trace.h_redsignal_raw;
    redsignal=redsignal/max(redsignal)*maxval;
    plot(xvals,redsignal,'Color',[153/256, 217/256, 234/256]);
    hold on;  plot(xvals,lefton,'k');
    hold on;  plot(xvals,righton,'m');    
    hold on;  plot(xvals,rep_data.trace.htrace,'Color',[0, 0, 1]);
    ylabel('H resp');

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

    nexttile(7,[1 2]);
    xvals=(1:length(rep_data.trace.vtrace))/10000;
    plot(xvals,rep_data.trace.vtrace);
    maxval=max(rep_data.trace.vtrace);
    downon=zeros(1,length(rep_data.trace.vtrace));
    downon(rep_data.trace.v_down_t)=maxval;
    upon=zeros(1,length(rep_data.trace.vtrace));
    upon(rep_data.trace.v_up_t)=maxval;
    redsignal=rep_data.trace.v_redsignal_raw;
    redsignal=redsignal/max(redsignal)*maxval;
    plot(xvals,redsignal,'Color',[153/256, 217/256, 234/256]);
    hold on;  plot(xvals,downon,'k');
    hold on;  plot(xvals,upon,'m');    
    hold on;  plot(xvals,rep_data.trace.vtrace,'Color',[0, 0, 1]);
    ylabel('V resp');

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
    elseif minRF==maxRF
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

    titlestr=[rep_data.strain,' ', rep_data.celltype,' ',rep_data.datestr, ' ', rep_data.cellid, ' ','rep_',num2str(rep_data.nrep_i)];
    title(t,titlestr, 'Interpreter', 'none');
    if ~isempty(resfolder)
%         if isempty(rep_data.cellid)
%             cellid='';
%         end
        if ~exist('extra_remark', 'var')
            figbasename=strjoin({'scanning_rep',rep_data.strain,rep_data.celltype, rep_data.datestr, rep_data.cellid,['rep',num2str(rep_data.nrep_i)]},'_');
        else
            figbasename=strjoin({'scanning_rep',rep_data.strain,rep_data.celltype, rep_data.datestr, rep_data.cellid,['rep',num2str(rep_data.nrep_i)], extra_remark},'_');
        end
        figbasename=[figbasename,'.fig'];
        figname = fullfile(resfolder,figbasename);
        savefig(figname);
    end
    close(f);
end