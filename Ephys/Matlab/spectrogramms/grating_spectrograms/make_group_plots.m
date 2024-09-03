% %make population plots
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';


%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
% group_labels=[["HS","shakB2"];["HSN","FlpND"];["HSN","FlpD"]]; 
group_labels=[["HS","FlpND"];["HS","FlpD"]]; 
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


file_ending='*grating_pd_nd_pwspectrum_*.mat';

%load all files that contain labels
%specified in the group
%get list of all files and folders 
path_patt=fullfile(population_folder,fullfile('**',file_ending));
filelist = dir(path_patt);

%% init the structure
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).pow_pd=[];
    mean_values_all(i).pow_nd=[];
    mean_values_all(i).pow_blpd_a=[];
    mean_values_all(i).pow_blpd_aa=[];
    mean_values_all(i).pow_blnd_a=[];
    mean_values_all(i).pow_blnd_aa=[];    
    mean_values_all(i).pwspectrum_pd=[];    
    mean_values_all(i).pwspectrum_nd=[];    
    mean_values_all(i).pwspectrum_blpd_a=[];    
    mean_values_all(i).pwspectrum_blpd_aa=[];    
    mean_values_all(i).pwspectrum_blnd_a=[];    
    mean_values_all(i).pwspectrum_blnd_aa=[];    
    mean_values_all(i).files={};
end
freq_vals=[];

%% iterate thru files 
% collect data only from the files that match filters and group patterns
for li=1:length(filelist)  
    filepath = filelist(li).folder;
    if contains(filepath,'compromised')
        continue;
    end
    
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
    if isempty(freq_vals)
        freq_vals= exp_data.pw.freq;
    end
        
    if isempty(mean_values_all(group_id).pwspectrum_pd)  
        i=1;
%    mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};    
    else
        i=size(mean_values_all(group_id).pwspectrum_pd,1)+1;
%         mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
    end    

    mean_values_all(group_id).pow_pd(i,:)=exp_data.pw.pow_pd_av;
    mean_values_all(group_id).pow_nd(i,:)=exp_data.pw.pow_nd_av;
    mean_values_all(group_id).pow_blpd_a(i,:)=exp_data.pw.pow_blpd_a_av;
    mean_values_all(group_id).pow_blpd_aa(i,:)=exp_data.pw.pow_blpd_aa_av;
    mean_values_all(group_id).pow_blnd_a(i,:)=exp_data.pw.pow_blnd_a_av;
    mean_values_all(group_id).pow_blnd_aa(i,:)=exp_data.pw.pow_blnd_aa_av;
   
    mean_values_all(group_id).pwspectrum_pd(i,:)=exp_data.pw.powerspectrum_pd_av;    
    mean_values_all(group_id).pwspectrum_nd(i,:)=exp_data.pw.powerspectrum_nd_av;    
    mean_values_all(group_id).pwspectrum_blpd_a(i,:)=exp_data.pw.powerspectrum_blpd_a_av; 
    mean_values_all(group_id).pwspectrum_blpd_aa(i,:)=exp_data.pw.powerspectrum_blpd_aa_av; 
    mean_values_all(group_id).pwspectrum_blnd_a(i,:)=exp_data.pw.powerspectrum_blnd_a_av; 
    mean_values_all(group_id).pwspectrum_blnd_aa(i,:)=exp_data.pw.powerspectrum_blnd_aa_av;     

    mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};    
end

npts=numel(freq_vals);
mean_spectrum_pd=zeros(ngrps,npts);
mean_spectrum_nd=zeros(ngrps,npts);
mean_spectrum_blpd_a=zeros(ngrps,npts);
mean_spectrum_blpd_aa=zeros(ngrps,npts);
mean_spectrum_blnd_a=zeros(ngrps,npts);
mean_spectrum_blnd_aa=zeros(ngrps,npts);
std_spectrum_pd=zeros(ngrps,npts);
std_spectrum_nd=zeros(ngrps,npts);
std_spectrum_blpd_a=zeros(ngrps,npts);
std_spectrum_blpd_aa=zeros(ngrps,npts);
std_spectrum_blnd_a=zeros(ngrps,npts);
std_spectrum_blnd_aa=zeros(ngrps,npts);

%for each group: low, mid, total
mean_pow_pd=zeros(ngrps,3);
std_pow_pd=zeros(ngrps,3);
mean_pow_nd=zeros(ngrps,3);
std_pow_nd=zeros(ngrps,3);
mean_pow_blpd_a=zeros(ngrps,3);
std_pow_blpd_a=zeros(ngrps,3);
mean_pow_blpd_aa=zeros(ngrps,3);
std_pow_blpd_aa=zeros(ngrps,3);
mean_pow_blnd_a=zeros(ngrps,3);
std_pow_blnd_a=zeros(ngrps,3);
mean_pow_blnd_aa=zeros(ngrps,3);
std_pow_blnd_aa=zeros(ngrps,3);

for i=1:ngrps
    mean_spectrum_pd(i,:)=mean(mean_values_all(i).pwspectrum_pd,1);
    std_spectrum_pd(i,:)=std(mean_values_all(i).pwspectrum_pd,[],1);
    mean_spectrum_nd(i,:)=mean(mean_values_all(i).pwspectrum_nd,1);
    std_spectrum_nd(i,:)=std(mean_values_all(i).pwspectrum_nd,[],1);
    mean_spectrum_blpd_a(i,:)=mean(mean_values_all(i).pwspectrum_blpd_a,1);
    std_spectrum_blpd_a(i,:)=std(mean_values_all(i).pwspectrum_blpd_a,[],1);
    mean_spectrum_blpd_aa(i,:)=mean(mean_values_all(i).pwspectrum_blpd_aa,1);
    std_spectrum_blpd_aa(i,:)=std(mean_values_all(i).pwspectrum_blpd_aa,[],1);
    mean_spectrum_blnd_a(i,:)=mean(mean_values_all(i).pwspectrum_blnd_a,1);
    std_spectrum_blnd_a(i,:)=std(mean_values_all(i).pwspectrum_blnd_a,[],1);
    mean_spectrum_blnd_aa(i,:)=mean(mean_values_all(i).pwspectrum_blnd_aa,1);
    std_spectrum_blnd_aa(i,:)=std(mean_values_all(i).pwspectrum_blnd_aa,[],1);

    mean_pow_pd(i,:)=mean(mean_values_all(i).pow_pd);
    std_pow_pd(i,:)=std(mean_values_all(i).pow_pd);
    mean_pow_nd(i,:)=mean(mean_values_all(i).pow_nd);
    std_pow_nd(i,:)=std(mean_values_all(i).pow_nd);
    mean_pow_blpd_a(i,:)=mean(mean_values_all(i).pow_blpd_a);
    std_pow_blpd_a(i,:)=std(mean_values_all(i).pow_blpd_a);
    mean_pow_blpd_aa(i,:)=mean(mean_values_all(i).pow_blpd_aa);
    std_pow_blpd_aa(i,:)=std(mean_values_all(i).pow_blpd_aa);
    mean_pow_blnd_a(i,:)=mean(mean_values_all(i).pow_blnd_a);
    std_pow_blnd_a(i,:)=std(mean_values_all(i).pow_blnd_a);
    mean_pow_blnd_aa(i,:)=mean(mean_values_all(i).pow_blnd_aa);
    std_pow_blnd_aa(i,:)=std(mean_values_all(i).pow_blnd_aa);     
end

legend_name={};
for i=1:ngrps    
    legend_name{i}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).pwspectrum_pd,1)),')'];
end

strain_colors = brewermap(ngrps,'Set1');

%% PD motion
figure, 
ymin=min([min(mean_spectrum_pd(:)),min(mean_spectrum_blpd_a(:)),min(mean_spectrum_blpd_aa(:))]);
ymax=max([max(mean_spectrum_pd(:)),max(mean_spectrum_blpd_a(:)),max(mean_spectrum_blpd_aa(:))]);

subplot(6,3,[1,2,3]);
for i=1:ngrps
    loglog(freq_vals, mean_spectrum_pd(i,:),"Color",strain_colors(i,:));
    hold on;
end
ylim([ymin,ymax]);
ylabel('PD motion');
legend(legend_name);

subplot(6,3,[4:6]);
for i=1:ngrps
    loglog(freq_vals, mean_spectrum_blpd_a(i,:),"Color",strain_colors(i,:));
    hold on;
end
ylim([ymin,ymax]);
ylabel('bl after PD');

subplot(6,3,[7:9]);
for i=1:ngrps
    loglog(freq_vals, mean_spectrum_blpd_aa(i,:),"Color",strain_colors(i,:));
    hold on;
end
ylim([ymin,ymax]);
ylabel('bl');

ytitlearr = {"< 10Hz", "[10,50] Hz", "all"};
for si=1:3
    subplot(6,3,10+(si-1)*3:10+si*3-1);

    b = bar([mean_pow_pd(:,si), mean_pow_blpd_a(:,si), mean_pow_blpd_aa(:,si)]','w');
    hold on;
    for i=1:ngrps
        nrec = size(mean_values_all(i).pwspectrum_pd,1);
        xdev =rand(1,nrec);
        xdev =repelem(xdev,1);
        barw=b(i).BarWidth*0.2;
        xdev=xdev*barw-barw/2;

        xvals=b(i).XEndPoints(1);
        xvals=repmat(xvals,[1,nrec]);
        xvals2plot=xvals+xdev;
        yvals = reshape(mean_values_all(i).pow_pd(:,si),1,[]);
        scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
        hold on;

        xvals=b(i).XEndPoints(2);
        xvals=repmat(xvals,[1,nrec]);
        xvals2plot=xvals+xdev;
        yvals = reshape(mean_values_all(i).pow_blpd_a(:,si),1,[]);
        scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
        hold on;

        xvals=b(i).XEndPoints(3);
        xvals=repmat(xvals,[1,nrec]);
        xvals2plot=xvals+xdev;
        yvals = reshape(mean_values_all(i).pow_blpd_aa(:,si),1,[]);
        scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
        hold on;
    end
    if si==3
        xlabel_name={"PD", "bl after PD", "bl"};
        xticklabels(xlabel_name);
        xtickangle(45);
    end
    ylabel(ytitlearr{si});
end


%% figure for ND direction
figure, 
ymin=min([min(mean_spectrum_nd(:)),min(mean_spectrum_blnd_a(:)),min(mean_spectrum_blnd_aa(:))]);
ymax=max([max(mean_spectrum_nd(:)),max(mean_spectrum_blnd_a(:)),max(mean_spectrum_blnd_aa(:))]);

subplot(6,3,[1,2,3]);
for i=1:ngrps
    loglog(freq_vals, mean_spectrum_nd(i,:),"Color",strain_colors(i,:));
    hold on;
end
ylim([ymin,ymax]);
ylabel('ND motion');
legend(legend_name);

subplot(6,3,[4:6]);
for i=1:ngrps
    loglog(freq_vals, mean_spectrum_blnd_a(i,:),"Color",strain_colors(i,:));
    hold on;
end
ylim([ymin,ymax]);
ylabel('bl after ND');

subplot(6,3,[7:9]);
for i=1:ngrps
    loglog(freq_vals, mean_spectrum_blnd_aa(i,:),"Color",strain_colors(i,:));
    hold on;
end
ylim([ymin,ymax]);
ylabel('bl');

ytitlearr = {"< 10Hz", "[10,50] Hz", "all"};
for si=1:3
    subplot(6,3,10+(si-1)*3:10+si*3-1);

    b = bar([mean_pow_nd(:,si), mean_pow_blnd_a(:,si), mean_pow_blnd_aa(:,si)]','w');
    hold on;
    for i=1:ngrps
        nrec = size(mean_values_all(i).pwspectrum_nd,1);
        xdev =rand(1,nrec);
        xdev =repelem(xdev,1);
        barw=b(i).BarWidth*0.2;
        xdev=xdev*barw-barw/2;

        xvals=b(i).XEndPoints(1);
        xvals=repmat(xvals,[1,nrec]);
        xvals2plot=xvals+xdev;
        yvals = reshape(mean_values_all(i).pow_nd(:,si),1,[]);
        scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
        hold on;

        xvals=b(i).XEndPoints(2);
        xvals=repmat(xvals,[1,nrec]);
        xvals2plot=xvals+xdev;
        yvals = reshape(mean_values_all(i).pow_blnd_a(:,si),1,[]);
        scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
        hold on;

        xvals=b(i).XEndPoints(3);
        xvals=repmat(xvals,[1,nrec]);
        xvals2plot=xvals+xdev;
        yvals = reshape(mean_values_all(i).pow_blnd_aa(:,si),1,[]);
        scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
        hold on;
    end
    if si==3
        xlabel_name={"ND", "bl after ND", "bl"};
        xticklabels(xlabel_name);
        xtickangle(45);
    end
    ylabel(ytitlearr{si});
end

