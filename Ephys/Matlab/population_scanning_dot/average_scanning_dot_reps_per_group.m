%go thru each folder in the parent folder and average all reps of scanning dots
% function average_scanning_dot_recording(parent_folder)

parent_folder = 'C:\Users\rsatapat\Documents\Victoria\RepositoryForScanningRect';
resfolder = 'C:\Users\rsatapat\Documents\Victoria\RepositoryForScanningRect\res';
% group_labels=["HSE","FlpD"];
% group_labels=["HSE","FlpND"];
group_labels=["HSE","FlpNDxDB331"];

% group_labels=[["HSN","FlpNDxDB331"]];%["HSE","FlpND DB331"]];
% group_labels=[["HSE","\FlpND\"]];
hor_degree_span_half=70;

%get all files with averages per cell
filelist=  get_all_scanning_dot_reps(parent_folder, group_labels);
nfiles = length(filelist);
group_data={};

%concatenate the data
nrec=length(filelist);
for li =1:nrec
    file_i=filelist{li};
    clear folder_data;

    %get one rep results
    load(file_i);
    group_data.cell(li)=folder_data;
end

%%resize all RF to the number of scan lines before averaging
[maxrow,maxcol] = get_RF_max_size(group_data);
[nrow, ncol] = get_RF_max_lines(group_data);
nrow=9;
ncol=13;
%get the grid corresponding to the scanlines
xvals=1:ncol;
yvals=1:nrow;
[xmesh,ymesh]=meshgrid(xvals,yvals);

%get [el,az] values corresponding to the grid of scanlines
ver_degree_span_half = 342*hor_degree_span_half/608;
azvals=linspace(-hor_degree_span_half,hor_degree_span_half,ncol);
elvals=linspace(-ver_degree_span_half,ver_degree_span_half,nrow);
[azmesh,elmesh] = meshgrid(azvals,elvals);

%wrap the grid on the sphere
%get [x,y] coordinatas in 3D corresponding to the grid of scanlines
azi=pi*azmesh/180;
eli=pi*elmesh/180;
[z3d,y3d,x3d]=sph2cart(eli,azi, 1);
spherical_coords_grid=cat(3,round(azi*180/pi),round(eli*180/pi));   
spherical_coords_grid_cart=cat(3,x3d,y3d,z3d);

%for each cell, downscale to the number of the scanlines,
% deform, then average
deformedRF = zeros(nrow,ncol,3,nrec);
nrec_real=0;
for ri=1:nrec
    RFxi=group_data.cell(ri).RF_x;
    if any(isnan(RFxi(:)))
        continue;
    end
    RFyi=group_data.cell(ri).RF_y;
    if any(isnan(RFyi(:)))
        continue;
    end
    nrec_real=nrec_real+1;

    %resize to max size
%     RFxi=imresize(RFxi,[nrow,ncol],"nearest");
%     RFyi=imresize(RFyi,[nrow,ncol],"nearest"); %% maybe use bicubic?
    RFxi=imresize(RFxi,[nrow,ncol],"bicubic");
    RFyi=imresize(RFyi,[nrow,ncol],"bicubic"); 

    
    %normalize, s.t. the max arrow is 1    
    lenmax=max((RFxi.^2+RFyi.^2).^0.5,[],"all");
    RFxi=RFxi./lenmax;
    RFyi=RFyi./lenmax;
    %deform onto a sphere
    
    %find the end of the arrow in 2D image
    xn=xmesh+RFxi;
    yn=ymesh+RFyi;
    %get corresponding az end el values
    azi=pi*interp1(xvals,azvals,xn,"linear","extrap")/180;
    eli=pi*interp1(yvals,elvals,yn,"linear","extrap")/180;   
    %get 3d coordinate of the end of the arrow
    [z3d,yn3d,xn3d]=sph2cart(eli,azi, 1);
    deformedRF(:,:,:,nrec_real) = cat(3,xn3d,yn3d,z3d)-spherical_coords_grid_cart;    
end

%average RF, angle and amplitude std
deformedRF(:,:,:,nrec_real+1:end)=[];
RF_x=squeeze(mean(deformedRF(:,:,1,:),4));
RF_y=squeeze(mean(deformedRF(:,:,2,:),4));
deformedRF_2D = deformedRF(:,:,1:2,:);
RF_amp = squeeze(vecnorm(deformedRF_2D,2,3));
RF_ampmean=mean(RF_amp,3);
RF_ampstd=std(RF_amp,[],3);
RF_alpha = squeeze(atan2(deformedRF_2D(:,:,1,:),deformedRF_2D(:,:,2,:)));
RF_alphamean=mean(unwrap(RF_alpha),3);
RF_alphastd=std(unwrap(RF_alpha),[],3);

%hor and vert vector profiles of the mean RF
RF_hor_vec_profile = cat(2,mean(RF_x,1)',mean(RF_y,1)');
RF_ver_vec_profile = cat(2,mean(RF_x,2),mean(RF_y,2));
RF_hor_profile = vecnorm(RF_hor_vec_profile,2,2);
RF_ver_profile = vecnorm(RF_ver_vec_profile,2,2);

%for each cell compute normalized horizontal profile
hor_profiles_all=zeros(nrec_real,length(RF_hor_profile));
hor_profiles_all_norm=zeros(nrec_real,length(RF_hor_profile));
for ci=1:nrec_real
    hprof=cat(2,mean(deformedRF(:,:,1,ci),1)',mean(deformedRF(:,:,2,ci),1)');
    hor_profiles_all(ci,:) = vecnorm(hprof,2,2);
    maxval=max(abs(hor_profiles_all(ci,:)));
    hor_profiles_all_norm(ci,:)=hor_profiles_all(ci,:)/maxval;   
end

%for each cell compute normalized vertical profile
ver_profiles_all=zeros(nrec_real,length(RF_ver_profile));
ver_profiles_all_norm=zeros(nrec_real,length(RF_ver_profile));
for ci=1:nrec_real
    vprof=cat(2,mean(deformedRF(:,:,1,ci),2),mean(deformedRF(:,:,2,ci),2));
    ver_profiles_all(ci,:) = vecnorm(vprof,2,2);
    maxval=max(abs(ver_profiles_all(ci,:)));
    ver_profiles_all_norm(ci,:)=ver_profiles_all(ci,:)/maxval;   
end

%for each cell compute the area of the contours
% some contours can be open.
% to close the the contours: padd the visual field by 1/100th 
% of min distance btw neiboring pts in the RF

%get the border and the area of the domain
dx=min(abs(diff(spherical_coords_grid_cart(1,:,1))))/100;
dy=min(abs(diff(spherical_coords_grid_cart(:,1,2))))/100;

paddded_xvals = padarray(spherical_coords_grid_cart(:,:,1),[1,1],0,'both');
paddded_xvals(:,1)=paddded_xvals(:,2)-dx;
paddded_xvals(:,end)=paddded_xvals(:,end-1)+dx;
paddded_xvals(1,:)=paddded_xvals(2,:);
paddded_xvals(end,:)=paddded_xvals(end-1,:);

paddded_yvals = padarray(spherical_coords_grid_cart(:,:,2),[1,1],0,'both');
paddded_yvals(1,:)=paddded_yvals(2,:)-dy;
paddded_yvals(end,:)=paddded_yvals(end-1,:)+dy;
paddded_yvals(:,end)=paddded_yvals(:,end-1);
paddded_yvals(:,1)=paddded_yvals(:,2);

xborder_p = [paddded_xvals(1,:), paddded_xvals(2:end-1,end)', paddded_xvals(end,end:-1:1), paddded_xvals(end-1:-1:1,1)']; 
yborder_p = [paddded_yvals(1,:), paddded_yvals(2:end-1,end)', paddded_yvals(end,end:-1:1), paddded_yvals(end-1:-1:1,1)']; 
border_area=polyarea(xborder_p,yborder_p);


cntr_areas=zeros(nrec_real,10-1);
cntr_frac_of_max=0.1:0.1:0.9;
cntrlevels=linspace(0, 0.9, 10);
for ci=1:nrec_real
    maxval=max(RF_amp(:,:,ci),[],"all");    
    RFamp_padded=RF_amp(:,:,ci)./maxval;
    RFamp_padded=padarray(RFamp_padded,[1,1],0,'both');
    %plotting RF for all cells
    f=figure; [mc, ~]= contourf(paddded_xvals,paddded_yvals,RFamp_padded,cntrlevels);
    [~,figname,~]=fileparts(filelist{ci});
    title(figname,'Interpreter','none');
    figname=fullfile(resfolder,[figname,'.pdf']);
    saveas(f,figname); close(f);
    %get contours without plotting
    % [mc, ~]= contour(paddded_xvals,paddded_yvals,RFamp_padded,cntrlevels);
    %get cumulative area of all superlevel sets
%     [~, level_cumarea] =get_total_contour_area(mc,xborder_p, yborder_p);
    [~, level_cumarea] =get_total_contour_area(mc);    
    cntr_areas(ci,:)=level_cumarea(:,2)/border_area;
end

%% find the areas of the superlevel sets on the mean RF
RF_amp_mean=(RF_x.^2+RF_y.^2).^0.5;
maxval=max(RF_amp_mean,[],"all");    
RFamp_padded=RF_amp_mean./maxval;
RFamp_padded=padarray(RFamp_padded,[1,1],0,'both');
f=figure; [mc, ~]= contourf(paddded_xvals,paddded_yvals,RFamp_padded,cntrlevels);
figname=strjoin(group_labels,'_');
figname=char(strrep(figname,'\',''));
title(figname,'Interpreter','none');
figname=fullfile(resfolder,[figname,'_avRF_cntr.pdf']);
saveas(f,figname); close(f);
    
% [~, level_cumarea] =get_total_contour_area(mc,xborder_p, yborder_p);
[~, level_cumarea] =get_total_contour_area(mc);
level_cumarea_meanRF=level_cumarea(:,2)/border_area;

%% save all in the data structure                
grp_data.grp_label = group_labels;
grp_data.files = filelist;
grp_data.RFall=deformedRF;
grp_data.deformed_grid=spherical_coords_grid_cart(:,:,1:2);
grp_data.spherical_grid=spherical_coords_grid;
grp_data.hspan_half=hor_degree_span_half;
grp_data.vspan_half=ver_degree_span_half;
grp_data.meanRF = cat(3,RF_x,RF_y);
grp_data.RF_hor_vec_profile= RF_hor_vec_profile;
grp_data.RF_ver_vec_profile= RF_ver_vec_profile;
grp_data.RF_hor_profile=RF_hor_profile;
grp_data.RF_ver_profile=RF_ver_profile;
%contours and hor profiles
grp_data.cntr_levels = cntr_frac_of_max;
grp_data.cntr_areas = cntr_areas;
grp_data.cntr_areas_meanRF = level_cumarea_meanRF;
grp_data.horizontal_profiles = hor_profiles_all;
grp_data.horizontal_profiles_normed1 = hor_profiles_all_norm;
grp_data.vertical_profiles = ver_profiles_all;
grp_data.vertical_profiles_normed1 = ver_profiles_all_norm;
%mean and std amplitude and angles
grp_data.RF_ampmean=RF_ampmean;
grp_data.RF_ampstd=RF_ampstd;
grp_data.RF_alphamean=RF_alphamean;
grp_data.RF_alphastd=RF_alphastd;

resfile=strjoin(grp_data.grp_label,'_');
resfile=strjoin(['scanning_rect_grp',resfile],'_');
resfilename=char(strrep(resfile,'\',''));
resfile=fullfile(resfolder,[resfilename,'.mat']);
save(resfile,'grp_data','-v7.3');

%make group plot
makeRFfigure_with_histograms_group_plot(grp_data, resfolder);
makeRFfigure_group_plot_all_cells(grp_data,resfolder);

xvals=linspace(-hor_degree_span_half,hor_degree_span_half,ncol);
figure, plot(xvals,fliplr(hor_profiles_all_norm)');
hold on; plot(xvals,fliplr(mean(hor_profiles_all_norm)),'LineWidth',2);
hold on; plot(xvals,fliplr(RF_hor_profile')/max(RF_hor_profile),'--','LineWidth',2);
title(resfilename,'Interpreter','none');

xvals=cntr_frac_of_max;
figure, plot(xvals,cntr_areas');
hold on; plot(xvals,mean(cntr_areas),'LineWidth',2);
hold on; plot(xvals,level_cumarea_meanRF,'--','LineWidth',2);
title(resfilename,'Interpreter','none');

% 
% end



function filelist=  get_all_scanning_dot_reps(parent_folder, group_labels)
    filelist ={};
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'scaning_rect_*_cellav.mat'];
    filelist_all = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    npat=size(group_labels,1);
    %go thru the list and collect all folders with '_scanning_rect_' files
    lii=1;
    for li =1:length(filelist_all)
        fullname=fullfile(filelist_all(li).folder,filelist_all(li).name);
        for patti=1:npat
            if str_contains_strarray(fullname,strcat(group_labels(patti,:), '_'))
                filelist{lii}=fullname;
                lii=lii+1;
                break;
            end
        end
    end    
end

function TF =  str_contains_strarray(strname,strarray)
    label_len=size(strarray,2);
    TF=1;
    for i=1:label_len
        TF=TF && contains(strname,strarray(i));
    end
end

function subfolders = get_all_scanning_dot_subfolders(parent_folder)
    subfolders ={};
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'scaning_rect_*_rep_*.mat'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with '_scanning_rect_' files    
    for li =1:length(filelist)          
        subfolders{li}=filelist(li).folder;                           
    end
    subfolders=unique(subfolders);
end

function scanning_rep_files = get_all_scanning_rep_files(parent_folder)
    scanning_rep_files ={};
        
    %get the list of files in all folders and subfolders
    allfiles=['scaning_rect_*_rep_*.mat'];
    filelist = dir(fullfile(parent_folder, allfiles));%get list of files and folders in any subfolder
    
    %go thru the list and collect all '_scanning_rect_' files       
    for li =1:length(filelist)          
        scanning_rep_files{li}=filelist(li).name;                           
    end    
end

function [maxrow,maxcol] = get_RF_max_size(group_data)
    nrec=numel(group_data.cell);
    maxrow=0;
    maxcol=0;
    for ri=1:nrec
        maxrow=max(maxrow,size(group_data.cell(ri).RF_y,1));
        maxcol=max(maxcol,size(group_data.cell(ri).RF_y,2));
    end
end

function [maxHlines, maxVlines] = get_RF_max_lines(group_data)
    nrec=numel(group_data.cell);
    maxHlines=0;
    maxVlines=0;
    for ri=1:nrec
        maxHlines=max(maxHlines,group_data.cell(ri).maxHlines);
        maxVlines=max(maxVlines,group_data.cell(ri).maxVlines);
    end
end

function [countrour_areas, level_cumarea] =get_total_contour_area(cntr_matrix)
n = 0;
i = 1;
contour_list={};
sz = size(cntr_matrix,2); %all pts
nn(1) = cntr_matrix(2,1); %pts of the 1st cntr
cx = cntr_matrix(1,2:nn(1)+1);
cy = cntr_matrix(2,2:nn(1)+1);
countrour_areas(i,1)= cntr_matrix(1,1); %level of the 1st cntr
countrour_areas(i,2) = polyarea(cx,cy);
contour_list(i).x=cx;
contour_list(i).y=cy;
while n+nn(i)+i < sz
    n = n + nn(i);
    i = i + 1;
    nn(i)=cntr_matrix(2,n+i); %pts of the i-th cntr
    cx = cntr_matrix(1,n+i+1:n+nn(i)+i);
    cy = cntr_matrix(2,n+i+1:n+nn(i)+i);
    countrour_areas(i,1)= cntr_matrix(1,n+i); %level of the i-th cntr
    countrour_areas(i,2) = polyarea(cx,cy);
    contour_list(i).x=cx;
    contour_list(i).y=cy;
end

ulevels=unique(countrour_areas(:,1));
level_cumarea=zeros(length(ulevels),2);
nlevels=length(ulevels);
% for i=1:nlevels
%     idx=find(countrour_areas(:,1)==ulevels(i));
%     figure,
%     for cii=1:length(idx)
%         plot(contour_list(idx(cii)).x,contour_list(idx(cii)).y);
%         hold on;
%     end
%     title(num2str(ulevels(i)));
% end

%check all contour of the same level for inclusion
for i=1:nlevels
    idx=find(countrour_areas(:,1)==ulevels(i));
    %get indices of contours of the level and their areas
    if length(idx)==1
        level_cumarea(i,:)=[ulevels(i),countrour_areas(idx,2)];
    else
        %sort contours by areas
        inds_area=[idx';countrour_areas(idx,2)']';
        inds_area=sortrows(inds_area,2,'descend');
        level_cumarea(i,:)=[ulevels(i),inds_area(1,2)];
        parent_x=contour_list(inds_area(1,1)).x;
        parent_y=contour_list(inds_area(1,1)).y;
        %for all smaller contour components check if they are inside a
        %bigger one
        for cii=2:length(idx)
            poly_x=contour_list(inds_area(cii,1)).x;
            poly_y=contour_list(inds_area(cii,1)).y;
            %if included
            if all(inpolygon(poly_x,poly_y,parent_x,parent_y))
                %substract area: it's a hole
                level_cumarea(i,2)=level_cumarea(i,2)-inds_area(cii,2);
            else %disjoint, add area, add poly into the parent list
                level_cumarea(i,2)=level_cumarea(i,2)+inds_area(cii,2);
                parent_x=[parent_x,poly_x];
                parent_y=[parent_y,poly_y];
            end
        end
    end
end
end

