function percentile_areas = makeRFfigure_with_histograms_group_plot(grp_data, resfolder)
    %horizontal and vertical components of the dynamic RF
    OF_hor= grp_data.meanRF(:,:,1);
    OF_ver= grp_data.meanRF(:,:,2);
    RFamp=(OF_hor.^2+OF_ver.^2).^0.5;
    maxval=max(RFamp(:));
    ncells=size(grp_data.RFall,4);
    
    %number of rows and columns in the low-res DRF
    [nrow, ncol]=size(OF_hor);
    %how many degrees does the half screen span horizontally (used for axis labels)
    half_az_span_deg = grp_data.hspan_half;
    half_el_span_deg = round(grp_data.vspan_half);
    
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
    [~,maxind] = find(grp_data.deformed_grid(:,:,2)==max(grp_data.deformed_grid(:,:,2),[],"all"));
    [maxcol,maxrow]=ind2sub(size(grp_data.deformed_grid(:,:,2)),maxind);
    yvals=grp_data.deformed_grid(:,maxcol,2);
%     nvb=length(grp_data.RF_hor_profile);    
    barh(yvals,grp_data.RF_ver_profile);
    xlim([0,1.2*max(grp_data.RF_ver_profile(:))]);
    xticks([0,min(grp_data.RF_ver_profile(:)),max(grp_data.RF_ver_profile(:))]);
    yticks([]);    
    ylabel('elevation');
       
    
    %subplot: OF
    ax3 = subplot('Position',pos3);        
    quiver(grp_data.deformed_grid(:,:,1),grp_data.deformed_grid(:,:,2),OF_hor,OF_ver, 0,'k','LineWidth',1);
    hold on;
    colcontours = colormap('winter');
    colcontours = colcontours([50,200],:);
    level_vals=[0.5*maxval,0.75*maxval];
    [MC, hContour] = contourf(grp_data.deformed_grid(:,:,1),grp_data.deformed_grid(:,:,2),RFamp,level_vals);
    drawnow;
    hFills = hContour.FacePrims;  % array of TriangleStrip objects
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    for idx = 1 : numel(hFills)
        hFills(idx).ColorData(:) = [colcontours(idx,:)*255, 100];
    end
    set(gca,'XDir','reverse');
    %get top 25 and 50 % area
    [contour_areas,contour_totareas]  = get_contour_areas(MC);
    percentile_areas=contour_totareas(:,2)';
    
    minx=min(grp_data.deformed_grid(:,:,1),[],"all");
    maxx=max(grp_data.deformed_grid(:,:,1),[],"all");
    miny=min(grp_data.deformed_grid(:,:,2),[],"all");
    maxy=max(grp_data.deformed_grid(:,:,2),[],"all");
    midx=mean([minx,maxx]);
    midy=mean([miny,maxy]);
    xlim([minx,maxx]);
    ylim([miny,maxy]);    
    axis equal;
    xticks([minx,midx,maxx]);
    yticks([miny,midy,maxy]);
    yticklabels([-half_el_span_deg,0,half_el_span_deg]);
    xlabels=[half_az_span_deg,0,-half_az_span_deg];
    xticklabels(xlabels);
    linkaxes([ax1,ax3],'y');
    
    %vertical bars
    ax2 = subplot('Position',pos2);     
    
    bar(grp_data.deformed_grid(1,:,1),grp_data.RF_hor_profile');
    ax2.YAxisLocation = 'right';
    xticks([]);    
    yticks([0,min(grp_data.RF_hor_profile(:)),max(grp_data.RF_hor_profile(:))]);
    ylim([0,1.2*max(grp_data.RF_ver_profile(:))]);
    set(gca,'XDir','reverse'); 
    xlabel('azimuth');
    linkaxes([ax3,ax2],'x');
    
    titlestr=strjoin(grp_data.grp_label,' ');
    titlestr=strjoin([titlestr, ['(',num2str(ncells),')']],' ');
    sgtitle(titlestr);  
    if ~isempty(resfolder)
        if ~exist(resfolder,'dir')
            mkdir(resfolder);
        end
        figbase=char(strjoin(grp_data.grp_label,'_'));
        fig_name=[figbase,'_Rf.fig'];
        fig_name=strrep(fig_name,'\','');
        fig_name=fullfile(resfolder,fig_name);
        savefig(fig_name);
        pdf_name=[fig_name(1:end-3),'pdf'];
        saveas(fig,pdf_name);
    end
end