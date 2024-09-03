% %make population plots
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';


%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
group_labels=[["HS","shakB2"];["HSN","FlpND"];["HSN","FlpD"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"]]; 
% group_labels=[["VS1_4","FlpD"];["VS1_4","FlpND"];["VS1_4","CantonS"]; ["VS1_4","FlpND DB331"]]; 
% group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","CantonS"]; ["HSE","FlpND DB331"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"];["HSN","CantonS"]; ["HSN","FlpND DB331"]]; 

ngrps= size(group_labels,1);

for ci=1:ngrps
    file_patt_i='';
    for cii=1:length(group_labels(ci,:))-1
        file_patt_i=strcat(file_patt_i,strcat(group_labels(ci,cii),'_'));        
    end
    file_patt_i = strcat(file_patt_i,group_labels(ci,end));
    file_patt{ci}=char(file_patt_i);
end


file_ending='*_pwspectrum_*.mat';

%load all files that contain labels
%specified in the group
%get list of all files and folders 
path_patt=fullfile(population_folder,fullfile('**',file_ending));
filelist = dir(path_patt);

%collect data only from the files that match filters and group patterns
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).lowfreq_frac=[];
    mean_values_all(i).midfreq_frac=[];
    mean_values_all(i).pwspectrum=[];    
    mean_values_all(i).files={};
end

midfreq = [];

for li=1:length(filelist)  
    filepath = filelist(li).folder;
    
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
    %exp_data is now available
    
    if isempty(midfreq)
        midfreq = find(exp_data.pw.freq<=50,1,'last');
        freq10 = find(exp_data.pw.freq<=10,1,'last');
    end

    %from powerspectrum of all rep, compute fraction per mid frequencies
    nrepi=size(exp_data.pw.powerspectrum_all,1);
    midfrac=[];
    for ri=1:nrepi
        totsum=sum(exp_data.pw.powerspectrum_all(ri,:));
        midfrac_i=sum(exp_data.pw.powerspectrum_all(ri,freq10:midfreq));
        midfrac=[midfrac,midfrac_i/totsum];
    end       

    if isempty(mean_values_all(group_id).pwspectrum)  
        mean_values_all(group_id).lowfreq_frac=mean(exp_data.pw.powfrac_lowfft);
        mean_values_all(group_id).midfreq_frac=mean(midfrac);        
        mean_values_all(group_id).pwspectrum=exp_data.pw.powerspectrum_av;
        mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};
    else
        i=size(mean_values_all(group_id).pwspectrum,1)+1;
        mean_values_all(group_id).lowfreq_frac(i,:)=mean(exp_data.pw.powfrac_lowfft);        
        mean_values_all(group_id).midfreq_frac(i,:)=mean(midfrac);        
        mean_values_all(group_id).pwspectrum(i,:)=exp_data.pw.powerspectrum_av;
        mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
    end
end

npts=size(mean_values_all(1).pwspectrum,2);
mean_spectrum=zeros(ngrps,npts);
std_spectrum=zeros(ngrps,npts);
mean_low_frac=zeros(ngrps,1);
mean_mid_frac=zeros(ngrps,1);
std_low_frac=zeros(ngrps,1);

for i=1:ngrps
    mean_spectrum(i,:)=mean(mean_values_all(i).pwspectrum,1);
    std_spectrum(i,:)=std(mean_values_all(i).pwspectrum,[],1);
    mean_low_frac(i)=mean(mean_values_all(i).lowfreq_frac);
    std_low_frac(i)=std(mean_values_all(i).lowfreq_frac);
    mean_mid_frac(i)=mean(mean_values_all(i).midfreq_frac);
end
legend_name={};
for i=1:ngrps    
    legend_name{i}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).pwspectrum,1)),')'];
end

freq = exp_data.pw.freq;
strain_colors = brewermap(ngrps,'Set1');
figure, 
subplot(2,2,[1,2]);
for i=1:ngrps
    loglog(freq, mean_spectrum(i,:));
    hold on;
end
legend(legend_name);
title('Power per Frequency')

subplot(2,2,3);
b = bar(mean_low_frac,'w');
hold on;

xlabel_name={};
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,1);
    xvals=b.XEndPoints(i);
    xvals=repmat(xvals,[1,nrec]);
    barw=b.BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;
    yvals = reshape(mean_values_all(i).lowfreq_frac,1,[]);
%     scatter(xvals,yvals,'filled','k');    
    scatter(xvals,yvals,7,'filled','k');    
    hold on;    
    xlabel_name{i}=[char(strjoin(group_labels(i,:),' '))];
end
xticklabels(xlabel_name);
xtickangle(45);
title('Power fraction by frequencies under 10Hz');


subplot(2,2,4);
b = bar(mean_mid_frac,'w');
hold on;
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,1);
    xvals=b.XEndPoints(i);
    xvals=repmat(xvals,[1,nrec]);
    barw=b.BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;
    yvals = reshape(mean_values_all(i).midfreq_frac,1,[]);
    scatter(xvals,yvals,7,'filled','k');    
    hold on;        
end
xticklabels(xlabel_name);
xtickangle(45);
title('Power fraction by frequencies 10-50Hz');



