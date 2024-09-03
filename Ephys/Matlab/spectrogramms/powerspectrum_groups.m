% %make population plots
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';


%% define groups within which you want to average power spectra

%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
group_labels=[["VS1_4","FlpD"];["VS1_4","FlpND"];["VS1_4","FlpND DB331"]]; 
% group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","CantonS"]; ["HSE","FlpND DB331"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"];["HSN","CantonS"]; ["HSN","FlpND DB331"]]; 


file_pw='*_pwspectrum.mat';
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
path_patt=fullfile(population_folder,fullfile('**',file_pw));
filelist = dir(path_patt);

%collect data only from the files that match filters and group patterns
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).pd_ps=[];
    mean_values_all(i).nd_ps=[];    
    mean_values_all(i).bl_ps=[];    
    mean_values_all(i).total_ps=[];    
    mean_values_all(i).files={};
end

%set of frequencies to re-sample individual spectra
fset=0:0.5:200;

for li=1:length(filelist)        
    filepath=fullfile(filelist(li).folder,filelist(li).name);
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
    
    load(filepath);
    
    %% bring power spectrum vals on the same frequence sample set   
    %interpolate pd power spectrum    
    xg=power_data.powerspectrum.pd.f(1:power_data.powerspectrum.pd.f_cut_ind);
    yg=power_data.powerspectrum.pd.p_normed(1:power_data.powerspectrum.pd.f_cut_ind);
    pd_pwvals=interp1(xg,yg,fset);
    %interpolate nd power spectrum    
    xg=power_data.powerspectrum.nd.f(1:power_data.powerspectrum.nd.f_cut_ind);
    yg=power_data.powerspectrum.nd.p_normed(1:power_data.powerspectrum.nd.f_cut_ind);
    nd_pwvals=interp1(xg,yg,fset);
    %interpolate baseline power spectrum   
    xg=power_data.powerspectrum.baseline.f(1:power_data.powerspectrum.baseline.f_cut_ind);
    yg=power_data.powerspectrum.baseline.p_normed(1:power_data.powerspectrum.baseline.f_cut_ind);
    bl_pwvals=interp1(xg,yg,fset);  
    
    %interpolate total power spectrum    
    xg=power_data.powerspectrum.total.f(1:power_data.powerspectrum.total.f_cut_ind);
    yg=power_data.powerspectrum.total.p_normed(1:power_data.powerspectrum.total.f_cut_ind);
    total_pwvals=interp1(xg,yg,fset);    
   
    if isempty(mean_values_all(group_id).pd_ps)  
        mean_values_all(group_id).pd_ps=pd_pwvals;
        mean_values_all(group_id).nd_ps=nd_pwvals;    
        mean_values_all(group_id).bl_ps=bl_pwvals;    
        mean_values_all(group_id).total_ps=total_pwvals;    
        mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};
    else
        i=size(mean_values_all(group_id).pd_ps,1)+1;        
        mean_values_all(group_id).pd_ps(i,:)=pd_pwvals;
        mean_values_all(group_id).nd_ps(i,:)=nd_pwvals;    
        mean_values_all(group_id).bl_ps(i,:)=bl_pwvals;    
        mean_values_all(group_id).total_ps(i,:)=total_pwvals;   
        mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
    end
end

%make a plot of average values
%find indices of non empty groups
nonemptygrps =[];
for i=1:ngrps
    if ~isempty(mean_values_all(i).pd_ps)
        nonemptygrps=[nonemptygrps,i];
    end
end
nunneg=length(nonemptygrps);
nfr=length(fset);
pd_mean=zeros(nunneg,nfr);
pd_std=zeros(nunneg,nfr);
nd_mean=zeros(nunneg,nfr);
nd_std=zeros(nunneg,nfr);
bl_mean=zeros(nunneg,nfr);
bl_std=zeros(nunneg,nfr);
total_mean=zeros(nunneg,nfr);
total_std=zeros(nunneg,nfr);

for ii=1:nunneg
    i=nonemptygrps(ii);
    pd_mean(ii,:)=mean(mean_values_all(i).pd_ps,1,'omitNaN');
    pd_std(ii,:)=std(mean_values_all(i).pd_ps,[],1,'omitNaN');
    
    nd_mean(ii,:)=mean(mean_values_all(i).nd_ps,1,'omitNaN');
    nd_std(ii,:)=std(mean_values_all(i).nd_ps,[],1,'omitNaN');
    
    bl_mean(ii,:)=mean(mean_values_all(i).bl_ps,1,'omitNaN');
    bl_std(ii,:)=std(mean_values_all(i).bl_ps,[],1,'omitNaN');
    
    total_mean(ii,:)=mean(mean_values_all(i).total_ps,1,'omitNaN');
    total_std(ii,:)=std(mean_values_all(i).total_ps,[],1,'omitNaN');
end
legend_name={};
for ii=1:nunneg
    i=nonemptygrps(ii);
    legend_name{ii}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).pd_ps,1)),')'];
end

strain_colors = brewermap(nunneg,'Set1');
figure, 
subplot(4,1,1);
for ii=1:nunneg    
    i=nonemptygrps(ii);
    stdshade_mean_std(total_mean(i,:),total_std(i,:),0.3,strain_colors(ii,:),fset);    
    hold on;
end
xlim([0,100]);
legend(legend_name,'Interpreter','non');
title('total trace');

subplot(4,1,2);
for ii=1:nunneg
    i=nonemptygrps(ii);
    stdshade_mean_std(pd_mean(i,:),pd_std(i,:),0.3,strain_colors(ii,:),fset);
    hold on;
end
xlim([0,100]);
title('PD trace');

subplot(4,1,3);
for ii=1:nunneg
    i=nonemptygrps(ii);
    stdshade_mean_std(nd_mean(i,:),nd_std(i,:),0.3,strain_colors(ii,:),fset);
    hold on;
end
xlim([0,100]);
title('ND trace');

subplot(4,1,4);
for ii=1:nunneg
    i=nonemptygrps(ii);
    stdshade_mean_std(bl_mean(i,:),bl_std(i,:),0.3,strain_colors(ii,:),fset);
    hold on;
end
xlim([0,100]);
title('baseline trace');




% save('VS14_spf01.mat','mean_values_all');
% %save figure
% filename =['grating_',save_name,file_ending];
% figname_fig=fullfile(population_folder,[filename,'_av_dir_resp.fig']);
% figname_png=fullfile(population_folder,[filename,'_av_dir_resp.png']);
% savefig(figname_fig);
% saveas(gcf,figname_png,'png');

