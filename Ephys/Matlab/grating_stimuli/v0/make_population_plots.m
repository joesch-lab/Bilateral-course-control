% %make population plots
% population_file = '\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository\HSS_all_recordings.mat';
% population_file = '\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository\HSN_all_recordings.mat';
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\grating_stimuli_res';
celltype=["HSE"]; %could be any set of cell types e.g.["HSE", "HSS"] or ["HS"]
speed = 5;
sp_freq=0.01;
file_patt='';
for ci=1:length(celltype)-1
    file_patt=strcat(file_patt,strcat(celltype(ci),'_'));
end
file_patt = strcat(file_patt,celltype(end));
file_patt=char(file_patt);

file_ending=['_dx',num2str(speed),'_spf',num2str(sp_freq*100)];

sampling_rate=10000;
baseline_fraction=0.25;

%load all file that start with 'grating_stimuli' and contain celltype
filelist=dir(population_folder);
filelist=filelist(~[filelist.isdir]);
i=1;
for li=1:length(filelist)
    if ~startsWith(filelist(li).name,'grating_stimuli_')
        continue;
    end
    if ~contains(filelist(li).name,celltype)
        continue;
    end
    matfile=fullfile(population_folder,filelist(li).name);
    load(matfile);
    
    %now the structure fly_av_vals is available            
    if celltype_av_vals.speed~=speed
        continue;
    end
    if celltype_av_vals.sp_freq~=sp_freq
        continue;
    end
    
    average_vals(i).av_resp = celltype_av_vals.av_resp;
    average_vals(i).std_resp = celltype_av_vals.std_resp;
    average_vals(i).av_bl = celltype_av_vals.av_bl;
    average_vals(i).std_bl = celltype_av_vals.std_bl;
    average_vals(i).nflies = celltype_av_vals.n_flies;
    average_vals(i).strain_type= celltype_av_vals.strain_type;    
    average_vals(i).dir_vector= celltype_av_vals.dir_vector;
    
    average_vals(i).cell_type= celltype;    
    average_vals(i).speed= speed;
    average_vals(i).sp_freq= sp_freq;
    i=i+1;
    clear celltype_av_vals;
end

dir_vector= average_vals(1).dir_vector;

%number of directions in the experiment
ndir=length(average_vals(1).dir_vector);
%number of strains 
nstrains=length(average_vals);

mindur_res=Inf;
mindur_bl=Inf;
%find min commong duration
for i=1:nstrains
    mindur_res=min(mindur_res,size(average_vals(i).av_resp,2));
    mindur_bl=min(mindur_res,size(average_vals(i).av_bl,2));
end

%crop resp and bl to the min duration
for i=1:nstrains
    average_vals(i).av_resp=average_vals(i).av_resp(:,1:mindur_res);
    average_vals(i).std_resp=average_vals(i).std_resp(:,1:mindur_res);
    average_vals(i).av_bl=average_vals(i).av_bl(:,1:mindur_bl);
    average_vals(i).std_bl=average_vals(i).std_bl(:,1:mindur_bl);    
end


%compute average responses
av_resp_per_direction = zeros(nstrains,ndir);
for i=1:nstrains
    av_resp_per_direction(i,:)=mean(average_vals(i).av_resp,2);
    std_resp_per_direction(i,:)=std(average_vals(i).av_resp,[],2);
    std_err(i,:)=std_resp_per_direction(i,:)/average_vals(i).nflies;
end
%find PD based on the activity of the first strain
[maxV, pdind] =max(av_resp_per_direction(1,:));
pdval=dir_vector(pdind);
if ndir==8
    ndind=pdind+4;
    ndind=1 + mod(ndind-1, ndir);
end
ndval=dir_vector(ndind);


% make_pd_nd_raw_plots
% strain_colors = distinguishable_colors(nstrains);
strain_colors = brewermap(nstrains,'Set1');
alpha=0.3;

figure,
bl_len=length(average_vals(1).av_bl(pdind,:));
bl_len_taken=round(baseline_fraction*bl_len);

stimstart=bl_len_taken+1;
x=1:(bl_len_taken+length(average_vals(1).av_resp(pdind,:)));
x=(x-stimstart)/sampling_rate;
title_str= ['PD and ND responses for ', file_patt,' cells'];
sgtitle(title_str,'FontSize',13, 'FontWeight','bold');
%pd responses
subplot(2,1,1);

legendstr={};
plotline_id=[];
minval=Inf;
maxval=-Inf;
for i=1:nstrains
    m=[average_vals(i).av_bl(pdind,end-bl_len_taken+1:end), average_vals(i).av_resp(pdind,:)];
    s=[average_vals(i).std_bl(pdind,end-bl_len_taken+1:end), average_vals(i).std_resp(pdind,:)];
    minval=min(minval,min(m-s));
    maxval=max(maxval,max(m+s));
    [lineOut, fillOut] = stdshade_mean_std(m,s,alpha,strain_colors(i,:),x);
    plotline_id=[plotline_id, lineOut];    
    lstri=[average_vals(i).strain_type, ' (',num2str(average_vals(i).nflies),')'];
    legendstr{i}=lstri;
    hold on;
end
plot([0,0],[minval,maxval],'--k');
legend(plotline_id, legendstr);
ytitle_str= ['PD responses (', num2str(pdval),' deg)'];
ylabel(ytitle_str);

%nd responses
subplot(2,1,2);

legendstr={};
plotline_id=[];
minval=Inf;
maxval=-Inf;
for i=1:nstrains
    m=[average_vals(i).av_bl(ndind,end-bl_len_taken+1:end), average_vals(i).av_resp(ndind,:)];
    s=[average_vals(i).std_bl(ndind,end-bl_len_taken+1:end), average_vals(i).std_resp(ndind,:)];
    minval=min(minval,min(m-s));
    maxval=max(maxval,max(m+s));
    [lineOut, fillOut] = stdshade_mean_std(m,s,alpha,strain_colors(i,:),x);
    plotline_id=[plotline_id, lineOut];    
    lstri=[average_vals(i).strain_type, ' (',num2str(average_vals(i).nflies),')'];
    legendstr{i}=lstri;
    hold on;
end

plot([0,0],[minval,maxval],'--k');
legend(plotline_id, legendstr);
ytitle_str= ['ND responses (', num2str(ndval),' deg)'];
ylabel(ytitle_str);

%save figure
filename =['grating_',file_patt,file_ending];
figname_fig=fullfile(population_folder,[filename,'_pdndtrace.fig']);
figname_png=fullfile(population_folder,[filename,'_pdndtrace.png']);
savefig(figname_fig);
saveas(gcf,figname_png,'png');



% make_average_all_dir_bar_plots();
figure,
b=bar(av_resp_per_direction');
barlgd = legend(legendstr);
barlgd.AutoUpdate='off';
legend()
hold on;
for i=1:nstrains
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*std_err(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end

xticklabels(dir_vector);
title_str= ['Directional responses for ', file_patt,' cells'];
title(title_str);
%save figure
filename =['grating_',file_patt,file_ending];
figname_fig=fullfile(population_folder,[filename,'_av_dir_resp_grating.fig']);
figname_png=fullfile(population_folder,[filename,'_av_dir_resp_grating.png']);
savefig(figname_fig);
saveas(gcf,figname_png,'png');
