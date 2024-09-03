% %make population plots
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';
resfolder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\scanning_rect_res';

%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
% group_labels=[["VS1_4","FlpD"];["VS1_4","FlpND"];["VS1_4","CantonS"]; ["VS1_4","FlpND DB331"]]; 
% group_labels=[["HSE","FlpND"];["HSN","FlpND"];["HSS","FlpND"]];
% group_labels=[["HSE","FlpD"];["HSE","FlpND"]];%["HSE","FlpND DB331"]];
% group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","CantonS"]; ["HSE","FlpND DB331"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"];["HSN","CantonS"]; ["HSN","FlpND DB331"]]; 
group_labels = [["H2","FlpD"];];

file_part='scaning_rect_*.mat';
ngrps= size(group_labels,1);

file_patt={};
for ci=1:ngrps
    file_patt_i='';
    for cii=1:length(group_labels(ci,:))-1
        file_patt_i=strcat(file_patt_i,strcat(group_labels(ci,cii),'_'));        
    end
    file_patt_i = strcat(file_patt_i,group_labels(ci,end));
    file_patt{ci}=char(file_patt_i);
end


%load all file that start with 'scaning_rect_' and contain labels
%specified in the group
%get list of all files and folders 
path_patt=fullfile(population_folder,fullfile('**',file_part));
filelist = dir(path_patt);

%collect data only from the files that match filters and group patterns
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);    
    mean_values_all(i).RFx=[];    
    mean_values_all(i).RFy=[];   
    mean_values_all(i).RFx_pos=[];
    mean_values_all(i).RFx_neg=[];
    mean_values_all(i).RFy_pos=[];
    mean_values_all(i).RFy_neg=[];        
    mean_values_all(i).files={};
    mean_values_all(i).n_flies=0;
end

for li=1:length(filelist)  
    disp(li);
    [filepath,~,~] = fileparts(filelist(li).folder);
    
    %split folder into words by slashes
    file_words=strsplit(filepath,{'\','/'});
    %find the match with the group: group names should match file_words
    group_id=-1;
    for gi=1:ngrps
        glen_i=size(group_labels(gi,:),2);
        match_i=0;
        for gli=1:glen_i
            if any(strcmpi(file_words,group_labels(gi,gli)))
                match_i=match_i+1;
            end
        end
        if match_i==glen_i
            group_id=gi;
            break;
        end
    end  
    if group_id==-1
%         disp([num2str(li),': no matching group');
        continue;
    end
   
    
    load(fullfile(filelist(li).folder,filelist(li).name));
    %normalize RF so that the len of the max vector is 1
    lenRF=(exp_data.RF.av_RFx.^2+exp_data.RF.av_RFy.^2).^0.5;
    maxlen=max(lenRF(:));
    maxlen_xpos=max(abs(exp_data.RF.av_RFx_pos(:)));
    maxlen_xneg=max(abs(exp_data.RF.av_RFx_neg(:)));
    maxlen_ypos=max(abs(exp_data.RF.av_RFy_pos(:)));
    maxlen_yneg=max(abs(exp_data.RF.av_RFy_neg(:)));
    if mean_values_all(group_id).n_flies==0          
        mean_values_all(group_id).RFx=exp_data.RF.av_RFx/maxlen;
        mean_values_all(group_id).RFy=exp_data.RF.av_RFy/maxlen;
        mean_values_all(group_id).RFx_pos=exp_data.RF.av_RFx_pos/maxlen_xpos;
        mean_values_all(group_id).RFx_neg=exp_data.RF.av_RFx_neg/maxlen_xneg;
        mean_values_all(group_id).RFy_pos=exp_data.RF.av_RFy_pos/maxlen_ypos;
        mean_values_all(group_id).RFy_neg=exp_data.RF.av_RFy_neg/maxlen_yneg;
        
        mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).n_flies=1;
    else        
        mean_values_all(group_id).n_flies=mean_values_all(group_id).n_flies+1; 
        
        mean_values_all(group_id).RFx=mean_values_all(group_id).RFx+exp_data.RF.av_RFx/maxlen;
        mean_values_all(group_id).RFy=mean_values_all(group_id).RFy+exp_data.RF.av_RFy/maxlen; 
        mean_values_all(group_id).RFx_pos=mean_values_all(group_id).RFx_pos + exp_data.RF.av_RFx_pos/maxlen_xpos;
        mean_values_all(group_id).RFx_neg=mean_values_all(group_id).RFx_neg + exp_data.RF.av_RFx_neg/maxlen_xneg;
        mean_values_all(group_id).RFy_pos=mean_values_all(group_id).RFy_pos + exp_data.RF.av_RFy_pos/maxlen_ypos;
        mean_values_all(group_id).RFy_neg=mean_values_all(group_id).RFy_neg + exp_data.RF.av_RFy_neg/maxlen_yneg;
        
        i= mean_values_all(group_id).n_flies;
        mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
    end
end

%mean across all flies
for i=1:ngrps    
    mean_values_all(i).RFx=mean_values_all(i).RFx/mean_values_all(i).n_flies;
    mean_values_all(i).RFy=mean_values_all(i).RFy/mean_values_all(i).n_flies;        
    mean_values_all(i).RFx_pos=mean_values_all(i).RFx_pos/mean_values_all(i).n_flies;   
    mean_values_all(i).RFx_neg=mean_values_all(i).RFx_neg/mean_values_all(i).n_flies;   
    mean_values_all(i).RFy_pos=mean_values_all(i).RFy_pos/mean_values_all(i).n_flies;   
    mean_values_all(i).RFy_neg=mean_values_all(i).RFy_neg/mean_values_all(i).n_flies;   
end


%% individual plot for a group (expand and deform)
for i=1:ngrps
    titlestr=[strrep(file_patt{i},'_',' '),' (',num2str(mean_values_all(i).n_flies),')'];
    figbase = fullfile(resfolder,file_patt{i}); 
    make_deformed_RF_figures(mean_values_all(i).RFx,mean_values_all(i).RFy, 15, titlestr,figbase);
    makeRF_indiv_compfigure(mean_values_all(i).RFx_pos, mean_values_all(i).RFx_neg, mean_values_all(i).RFy_pos, mean_values_all(i).RFy_neg, 15, titlestr, figbase);
end

%% plot the groups together in one quiver plot
%create one figure with to plots: Optic flow 
fig=figure('Position',[100,100,800,600]); 
nbin=15;
%how many degrees does the half screen span horizontally (used for axis labels)
half_az_span_deg =70;

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

ax1 = subplot('Position',pos1);  
ax3 = subplot('Position',pos3);     
ax2 = subplot('Position',pos2);     
minbarvert=Inf;
maxbarvert=-Inf;
minbarhor=Inf;
maxbarhor=-Inf;
cmap=brewermap(ngrps,'set1');
maxlen=0;
%%for the clarity it isn't recommended to plot more than 2groups
grpoups_to_plot=[1,2];
maxleni=zeros(1,length(grpoups_to_plot));
for gi=1:length(grpoups_to_plot)
    i=grpoups_to_plot(gi);
    %make the largest arrow have length 1;
    [bin_dx_av, bin_resp, bin_resp_deformed, ...
        bin_hor, bin_ver, bin_hor_sep, bin_ver_sep, el_bincenters, az_bincenters] ...
        = deform_plot_uniform_sampling(mean_values_all(i).RFx,mean_values_all(i).RFy, 15, 'c', 0);
    
    %horizontal and vertical components of the dynamic RF
    OF_hor= bin_dx_av(:,:,1);
    OF_ver= bin_dx_av(:,:,2);

    leni=(OF_hor.^2+OF_ver).^0.5;
    maxleni(gi)=max(maxleni(gi),max(leni(:)));
end
for i=1:2%ngrps
    %deform RF
    [bin_dx_av, bin_resp, bin_resp_deformed, ...
        bin_hor, bin_ver, bin_hor_sep, bin_ver_sep, el_bincenters, az_bincenters] ...
        = deform_plot_uniform_sampling(mean_values_all(i).RFx/maxleni(gi),mean_values_all(i).RFy/maxleni(gi), 15, 'c', 0);
    
    %horizontal and vertical components of the dynamic RF
    OF_hor= bin_dx_av(:,:,1);
    OF_ver= bin_dx_av(:,:,2);
    
    %number of rows and columns in the low-res DRF
    [nrow, ncol,xy]=size(bin_dx_av);
    
    %horizontal bars    
    nvb=length(bin_ver);
    axes(ax1);
    hold on;
    barh(1:nvb,bin_ver', 'FaceColor',cmap(i,:),'FaceAlpha',0.3);
    minbarvert=min(minbarvert,min(bin_ver));
    maxbarvert=max(maxbarvert,max(bin_ver));
 
    
    %subplot: OF    
%     lengths=sqrt(OF_hor.^2.+OF_ver.^2);
%     lenmax=max(lengths(:));
%     OF_hor_plot=OF_hor./lenmax;
%     OF_ver_plot=OF_ver./lenmax;
    [x,y] = meshgrid(1:ncol, 1:nrow);
    axes(ax3);
    hold on;
    quiver(ax3,x,y,OF_hor,OF_ver, 0,'color', cmap(i,:),'LineWidth',1);         
    
    
    %vertical bars    
    nhb=length(bin_hor);
    axes(ax2);
    hold on;
    bar(1:nhb,bin_hor','FaceColor',cmap(i,:), 'FaceAlpha',0.3); 
    minbarhor=min(minbarhor,min(bin_hor));
    maxbarhor=max(maxbarhor,max(bin_hor));   
end

%labels, orientation, limits of axis

axes(ax1);
xlim([minbarvert,maxbarvert]);
xticks([minbarvert,maxbarvert]);
yticks([]);    
ylabel('elevation');


axes(ax3);
hold on;
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

axes(ax2);
hold on;
ax2.YAxisLocation = 'right';
xticks([]);    
ylim([minbarhor,maxbarhor]);
yticks([minbarhor,maxbarhor]);
set(gca,'XDir','reverse'); 
xlabel('azimuth');
linkaxes([ax3,ax2],'x');

name1=[strrep(file_patt{grpoups_to_plot(1)},'_',' '),' (',num2str(mean_values_all(grpoups_to_plot(1)).n_flies),')'];
name2=[strrep(file_patt{grpoups_to_plot(2)},'_',' '),' (',num2str(mean_values_all(grpoups_to_plot(2)).n_flies),')'];
sgtitle( sprintf('%s{%f %f %f}%s %s{%f %f %f}vs %s{%f %f %f}%s', '\color[rgb]',cmap(1,:), name1, '\color[rgb]',[0,0,0],  '\color[rgb]', cmap(2,:), name2), 'interpreter', 'tex');

figfile=[file_patt{grpoups_to_plot(1)},'_vs_',file_patt{grpoups_to_plot(2)},'.fig'];
figfile=fullfile(resfolder, figfile);
savefig(figfile);

