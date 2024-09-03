function makeRFfigure_group_plot_all_cells(grp_data, resfolder)
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
    fig = figure;    
    for ri=1:ncells
        quiver(grp_data.deformed_grid(:,:,1),grp_data.deformed_grid(:,:,2),grp_data.RFall(:,:,1,ri),grp_data.RFall(:,:,2,ri), 0,'k','LineWidth',0.25);
        hold on;
    end
    quiver(grp_data.deformed_grid(:,:,1),grp_data.deformed_grid(:,:,2),OF_hor,OF_ver, 0,'k','LineWidth',2);
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
    set(gca,'XDir','reverse');
    
    titlestr=strjoin(grp_data.grp_label,' ');
    titlestr=strjoin([titlestr, ['(',num2str(ncells),')']],' ');
    sgtitle(titlestr);  

    if ~isempty(resfolder)
        if ~exist(resfolder,'dir')
            mkdir(resfolder);
        end
        figbase=char(strjoin(grp_data.grp_label,'_'));
        figbase=strrep(figbase,'\','');
        fig_name=[figbase,'_Rf_allcells.fig'];        
        fig_name=fullfile(resfolder,fig_name);
        savefig(fig_name);
        pdf_name=[fig_name(1:end-3),'pdf'];
        saveas(fig,pdf_name);
    end
end
    