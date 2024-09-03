% %make population plots
% population_file = '\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository\HSS_all_recordings.mat';
% population_file = '\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository\HSN_all_recordings.mat';
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';
save_name = 'test_norm';
%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","Canton"]]; 
speed = 5;
sp_freq=0.01;
% option to normalize max response to of each group to 1
normalize_to_one = 0;

file_ending=['_dx',num2str(speed),'_spf',num2str(sp_freq*100)];

sampling_rate=10000;
baseline_fraction=0.25;

ngrps= size(group_labels,1);

for ci=1:ngrps
    file_patt_i='';
    for cii=1:length(group_labels(ci,:))-1
        file_patt_i=strcat(file_patt_i,strcat(group_labels(ci,cii),'_'));        
    end
    file_patt_i = strcat(file_patt_i,group_labels(ci,end));
    file_patt{ci}=char(file_patt_i);
end


%load all file that start with 'grating_stimuli' and contain labels
%specified in the group
%get list of all files and folders 
path_patt=fullfile(population_folder,fullfile('**','*.mat'));
filelist = dir(path_patt);

%remove files that do not contain grating_stimuli_res in the path and that
%do not start with the date
exclude_inds=[];

for li=1:length(filelist)
    if ~startsWith(filelist(li).name,digitsPattern(6))
        exclude_inds=[exclude_inds,li];
        continue;
    end
    if ~contains(filelist(li).folder,'grating_stimuli_res')
        exclude_inds=[exclude_inds,li];
        continue;
    end
end
filelist(exclude_inds)=[];

%for each defined group select files that mathch label and pool the average
%values
group_vals={};
mindur_resp=Inf;
mindur_bl=Inf;

for gi=1:ngrps
    group_vals(gi).group_name=file_patt{gi};
    
    %select files
    f_ids = get_folder_indices_with_pattern({filelist(:).folder},group_labels(gi,:));    
    i=1; 
    average_vals={};
    for fi=1:length(f_ids)
        matfile=fullfile(filelist(f_ids(fi)).folder,filelist(f_ids(fi)).name);
        load(matfile);

        %now the structure fly_av_vals is available            
        if fly_av_vals.speed~=speed
            continue;
        end
        if fly_av_vals.sp_freq~=sp_freq
            continue;
        end

        average_vals(i).av_resp = fly_av_vals.av_resp;
        average_vals(i).std_resp = fly_av_vals.std_resp;
        average_vals(i).av_bl = fly_av_vals.av_bl;
        average_vals(i).std_bl = fly_av_vals.std_bl;        
        average_vals(i).strain_type= fly_av_vals.strain_type;    
        average_vals(i).dir_vector= fly_av_vals.dir_vector;
        average_vals(i).cell_type= fly_av_vals.cell_type;    
        average_vals(i).speed= speed;
        average_vals(i).sp_freq= sp_freq;
        
        mindur_resp=min(mindur_resp,size(average_vals(i).av_resp,2));
        mindur_bl=min(mindur_bl,size(average_vals(i).av_bl,2));
        i=i+1;
        clear fly_av_vals;
    end
    group_vals(gi).all_vals=average_vals;
end

dir_vector= group_vals(1).all_vals(1).dir_vector;
%number of directions in the experiment
ndir=length(average_vals(1).dir_vector);

for gi=1:ngrps
    nflies = length(group_vals(gi).all_vals);
    all_flies_resp=zeros(nflies,ndir,mindur_resp);
    all_flies_bl=zeros(nflies,ndir,mindur_bl);
    for fi=1:nflies
        all_flies_resp(fi,:,:)=group_vals(gi).all_vals(fi).av_resp(:,1:mindur_resp);
        all_flies_bl(fi,:,:)=group_vals(gi).all_vals(fi).av_bl(:,1:mindur_bl);
    end
    group_vals(gi).nflies=nflies;
    group_vals(gi).av_resp=squeeze(mean(all_flies_resp,1));
    group_vals(gi).std_resp=squeeze(std(all_flies_resp,[],1));
    group_vals(gi).av_bl=squeeze(mean(all_flies_bl,1));
    group_vals(gi).std_bl=squeeze(std(all_flies_bl,[],1));
    %average the whole response
    group_vals(gi).av_per_dir=mean(all_flies_resp,[1,3]);
    group_vals(gi).std_av_per_dir=std(all_flies_resp,[],[1,3]);
    group_vals(gi).std_err=group_vals(gi).std_av_per_dir/nflies;
end


%find PD based on the activity of the first group
[maxV, pdind] =max(group_vals(1).av_per_dir);
pdval=dir_vector(pdind);
if ndir==8
    ndind=pdind+4;
    ndind=1 + mod(ndind-1, ndir);
end
ndval=dir_vector(ndind);


% make_pd_nd_raw_plots
strain_colors = brewermap(ngrps,'Set1');
alpha=0.3;

figure,
bl_len_taken=round(baseline_fraction*mindur_bl);

stimstart=bl_len_taken+1;
x=1:(bl_len_taken+mindur_resp);
x=(x-stimstart)/sampling_rate;
title_str= ['PD and ND responses'];
sgtitle(title_str,'FontSize',13, 'FontWeight','bold');
%pd responses
subplot(2,1,1);

legendstr={};
plotline_id=[];
minval=Inf;
maxval=-Inf;
for gi=1:ngrps
    m=[group_vals(gi).av_bl(pdind,end-bl_len_taken+1:end), group_vals(gi).av_resp(pdind,:)];
    s=[group_vals(gi).std_bl(pdind,end-bl_len_taken+1:end), group_vals(gi).std_resp(pdind,:)];
    minval=min(minval,min(m-s));
    maxval=max(maxval,max(m+s));
    if normalize_to_one
        max_val_normed=max(abs(m));
        m=m/max_val_normed;
        s=s/max_val_normed;           
    end
    [lineOut, fillOut] = stdshade_mean_std(m,s,alpha,strain_colors(gi,:),x);
    plotline_id=[plotline_id, lineOut];    
    lstri=[group_vals(gi).group_name, ' (',num2str(group_vals(gi).nflies),')'];
    legendstr{gi}=lstri;
    hold on;
end
if normalize_to_one        
    minval=minval/max_val_normed;    
    maxval=maxval/max_val_normed;        
end
plot([0,0],[minval,maxval],'--k');
legend(plotline_id, legendstr, 'Interpreter', 'none');
ytitle_str= ['PD responses (', num2str(pdval),' deg)'];
ylabel(ytitle_str);

%nd responses
subplot(2,1,2);

legendstr={};
plotline_id=[];
minval=Inf;
maxval=-Inf;
for gi=1:ngrps
    m=[group_vals(gi).av_bl(ndind,end-bl_len_taken+1:end), group_vals(gi).av_resp(ndind,:)];
    s=[group_vals(gi).std_bl(ndind,end-bl_len_taken+1:end), group_vals(gi).std_resp(ndind,:)];
    minval=min(minval,min(m-s));
    maxval=max(maxval,max(m+s));
    if normalize_to_one
        max_val_normed=max(abs(m));
        m=m/max_val_normed;
        s=s/max_val_normed;           
    end
    [lineOut, fillOut] = stdshade_mean_std(m,s,alpha,strain_colors(gi,:),x);
    plotline_id=[plotline_id, lineOut];    
    lstri=[group_vals(gi).group_name, ' (',num2str(group_vals(gi).nflies),')'];
    legendstr{gi}=lstri;
    hold on;
end
if normalize_to_one        
    minval=minval/max_val_normed;    
    maxval=maxval/max_val_normed;        
end
plot([0,0],[minval,maxval],'--k');
legend(plotline_id, legendstr, 'Interpreter', 'none');
ytitle_str= ['ND responses (', num2str(ndval),' deg)'];
ylabel(ytitle_str);

%save figure
filename =['grating_',save_name,file_ending];
figname_fig=fullfile(population_folder,[filename,'_pdndtrace.fig']);
figname_png=fullfile(population_folder,[filename,'_pdndtrace.png']);
savefig(figname_fig);
saveas(gcf,figname_png,'png');



% make_average_all_dir_bar_plots();
%place all averages per direction into one array
all_av_dir=zeros(ngrps,ndir);
for gi=1:ngrps
    all_av_dir(gi,:)=group_vals(gi).av_per_dir;
end

figure,
b=bar(all_av_dir');
barlgd = legend(legendstr, 'Interpreter', 'none');
barlgd.AutoUpdate='off';
legend()
hold on;
for i=1:ngrps
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*group_vals(i).std_err(:)';    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end

xticklabels(dir_vector);
title_str= ['Directional responses'];
title(title_str);
%save figure
filename =['grating_',save_name,file_ending];
figname_fig=fullfile(population_folder,[filename,'_av_dir_resp.fig']);
figname_png=fullfile(population_folder,[filename,'_av_dir_resp.png']);
savefig(figname_fig);
saveas(gcf,figname_png,'png');


function ids = get_folder_indices_with_pattern(folderlist,pattarr)
    np=size(pattarr,2);
    loginds=ones(1, size(folderlist,2));
    for pi=1:np
        yn=contains(folderlist,pattarr(pi));
        loginds = loginds& yn;
    end
    ids= find(loginds);
end