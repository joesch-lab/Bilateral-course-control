% %make population plots
clear all;
population_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';


%can define any group to plot, the labels of the groups
% are defined within [] and separated by ;
% group_labels=[["HS","shakB2"];["HSN","FlpND"];["HSN","FlpD"]];
group_labels=[["HS","FlpND"];["HS","FlpD"];["HS","shakB2"]];
% group_labels=[["VS","FlpND"];["VS","shakB2"]];
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


file_ending='*flashes_on_off_pwspectrum_*.mat';

%load all files that contain labels
%specified in the group
%get list of all files and folders
path_patt=fullfile(population_folder,fullfile('**',file_ending));
filelist = dir(path_patt);

%% init the structure
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).pow_on=[];
    mean_values_all(i).pow_off=[];
    mean_values_all(i).pwspectrum_on=[];
    mean_values_all(i).pwspectrum_off=[];
    mean_values_all(i).pow_on_normed=[];
    mean_values_all(i).pow_off_normed=[];
    mean_values_all(i).pwspectrum_on_normed=[];
    mean_values_all(i).pwspectrum_off_normed=[];
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
        f10=find(freq_vals<=10,1,'last');
        f50=find(freq_vals<=50,1,'last');
    end

    if isempty(mean_values_all(group_id).pwspectrum_on)
        i=1;
    else
        i=size(mean_values_all(group_id).pwspectrum_on,1)+1;
    end

    mean_values_all(group_id).pow_on(i,:)=exp_data.pw.pow_on_av;
    mean_values_all(group_id).pow_off(i,:)=exp_data.pw.pow_off_av;
    mean_values_all(group_id).pwspectrum_on(i,:)=exp_data.pw.powerspectrum_on_av;
    mean_values_all(group_id).pwspectrum_off(i,:)=exp_data.pw.powerspectrum_off_av;
    %normalized spectrum, s.t. total power is one
    mean_values_all(group_id).pwspectrum_on_normed(i,:)=exp_data.pw.powerspectrum_on_av/sum(exp_data.pw.powerspectrum_on_av);
    mean_values_all(group_id).pwspectrum_off_normed(i,:)=exp_data.pw.powerspectrum_off_av/sum(exp_data.pw.powerspectrum_off_av);
    mean_values_all(group_id).pow_on_normed(i,1)=sum(mean_values_all(group_id).pwspectrum_on_normed(i,1:f10));
    mean_values_all(group_id).pow_on_normed(i,2)=sum(mean_values_all(group_id).pwspectrum_on_normed(i,f10+1:f50));
    mean_values_all(group_id).pow_on_normed(i,3)=sum(mean_values_all(group_id).pwspectrum_on_normed(i,f50+1:end));
    mean_values_all(group_id).pow_off_normed(i,1)=sum(mean_values_all(group_id).pwspectrum_off_normed(i,1:f10));
    mean_values_all(group_id).pow_off_normed(i,2)=sum(mean_values_all(group_id).pwspectrum_off_normed(i,f10+1:f50));
    mean_values_all(group_id).pow_off_normed(i,3)=sum(mean_values_all(group_id).pwspectrum_off_normed(i,f50+1:end));

    mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
end

npts=numel(freq_vals);
mean_spectrum_on=zeros(ngrps,npts);
mean_spectrum_off=zeros(ngrps,npts);
std_spectrum_on=zeros(ngrps,npts);
std_spectrum_off=zeros(ngrps,npts);

%for each group: low, mid, total
mean_pow_on=zeros(ngrps,3);
std_pow_on=zeros(ngrps,3);
mean_pow_off=zeros(ngrps,3);
std_pow_off=zeros(ngrps,3);

%mean of normalized spectra
mean_spectrum_on_normed=zeros(ngrps,npts);
mean_spectrum_off_normed=zeros(ngrps,npts);
std_spectrum_on_normed=zeros(ngrps,npts);
std_spectrum_off_normed=zeros(ngrps,npts);
mean_pow_on_normed=zeros(ngrps,3);
std_pow_on_normed=zeros(ngrps,3);
mean_pow_off_normed=zeros(ngrps,3);
std_pow_off_normed=zeros(ngrps,3);



for i=1:ngrps
    mean_spectrum_on(i,:)=mean(mean_values_all(i).pwspectrum_on,1);
    std_spectrum_on(i,:)=std(mean_values_all(i).pwspectrum_on,[],1);
    mean_spectrum_off(i,:)=mean(mean_values_all(i).pwspectrum_off,1);
    std_spectrum_off(i,:)=std(mean_values_all(i).pwspectrum_off,[],1);

    mean_pow_on(i,:)=mean(mean_values_all(i).pow_on);
    std_pow_on(i,:)=std(mean_values_all(i).pow_on);
    mean_pow_off(i,:)=mean(mean_values_all(i).pow_off);
    std_pow_off(i,:)=std(mean_values_all(i).pow_off);

    mean_spectrum_on_normed(i,:)=mean(mean_values_all(i).pwspectrum_on_normed,1);
    std_spectrum_on_normed(i,:)=std(mean_values_all(i).pwspectrum_on_normed,[],1);
    mean_spectrum_off_normed(i,:)=mean(mean_values_all(i).pwspectrum_off_normed,1);
    std_spectrum_off_normed(i,:)=std(mean_values_all(i).pwspectrum_off_normed,[],1);

    mean_pow_on_normed(i,:)=mean(mean_values_all(i).pow_on_normed);
    std_pow_on_normed(i,:)=std(mean_values_all(i).pow_on_normed);
    mean_pow_off_normed(i,:)=mean(mean_values_all(i).pow_off_normed);
    std_pow_off_normed(i,:)=std(mean_values_all(i).pow_off_normed);
end

legend_name={};
for i=1:ngrps
    legend_name{i}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).pwspectrum_on,1)),')'];
end

strain_colors = brewermap(ngrps,'Set1');

%% make plots
% mean power decompostion
figure, t= tiledlayout(3,2,"TileSpacing","compact");
f_s = find(freq_vals>=1,1,'first');
f_e = find(freq_vals<=100,1,'last');
ymin=min([min(mean_spectrum_on(:,f_s:f_e),[],"all"),min(mean_spectrum_off(:,f_s:f_e),[],"all")]);
ymax=max([max(mean_spectrum_on(:,f_s:f_e),[],"all"),max(mean_spectrum_off(:,f_s:f_e),[],"all")]);

nexttile(1,[1,2]);%subplot(3,4,1:4);
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_on,1);    
    for ri=1:nrec        
        loglog(freq_vals(f_s:f_e), mean_values_all(i).pwspectrum_on(ri,f_s:f_e),"Color",[strain_colors(i,:),0.2]);
        hold on;
    end   
end
axa=[];
for i=1:ngrps    
    axai = loglog(freq_vals(f_s:f_e), mean_spectrum_on(i,f_s:f_e),"Color",strain_colors(i,:),'LineWidth',2);
    hold on;
    axa=[axa,axai];
end

ylim([ymin,ymax]);
ylabel('ON resp');
legend(axa, legend_name);

nexttile(3,[1,2]);%subplot(3,4,[5:8]);
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_on,1);    
    for ri=1:nrec        
        loglog(freq_vals(f_s:f_e), mean_values_all(i).pwspectrum_off(ri,f_s:f_e),"Color",[strain_colors(i,:),0.2]);
        hold on;
    end
    loglog(freq_vals(f_s:f_e), mean_spectrum_off(i,f_s:f_e),"Color",strain_colors(i,:),'LineWidth',2);
    hold on;
end
ylim([ymin,ymax]);
ylabel('OFF resp');

%bars of power contribution within each band
xlabelarr = {"< 10Hz", "[10,50] Hz", "all"};

nexttile;%subplot(3,4,[9:10]);
b = bar(mean_pow_on','w');
hold on;
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_on,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,1);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;

    xvals=b(i).XEndPoints(1);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_on(:,1),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(2);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_on(:,2),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(3);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_on(:,3),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;
end
xticklabels(xlabelarr );
xtickangle(45);
ylabel('ON');


nexttile;%subplot(3,4,[11:12]);
b = bar(mean_pow_off','w');
hold on;
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_off,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,1);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;

    xvals=b(i).XEndPoints(1);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_off(:,1),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(2);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_off(:,2),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(3);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_off(:,3),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;
end
xticklabels(xlabelarr );
xtickangle(45);
ylabel('OFF');
title(t,'raw power spectrum');



%% make plots - normalized values
% mean power decompostion
figure, t= tiledlayout(3,2,"TileSpacing","compact");
ymin=min([min(mean_spectrum_on_normed(:,f_s:f_e),[],"all"),min(mean_spectrum_off_normed(:,f_s:f_e),[],"all")]);
ymax=max([max(mean_spectrum_on_normed(:,f_s:f_e),[],"all"),max(mean_spectrum_off_normed(:,f_s:f_e),[],"all")]);

nexttile(1,[1,2]);%subplot(3,4,1:4);
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_on_normed,1);    
    for ri=1:nrec        
        loglog(freq_vals(f_s:f_e), mean_values_all(i).pwspectrum_on_normed(ri,f_s:f_e),"Color",[strain_colors(i,:),0.2]);
        hold on;
    end   
end
axa=[];
for i=1:ngrps    
    axai = loglog(freq_vals(f_s:f_e), mean_spectrum_on_normed(i,f_s:f_e),"Color",strain_colors(i,:),'LineWidth',2);
    hold on;
    axa=[axa,axai];
end
ylim([ymin,ymax]);
ylabel('ON resp');
legend(axa, legend_name);

nexttile(3,[1,2]);%subplot(3,4,[5:8]);
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_on_normed,1);    
    for ri=1:nrec        
        loglog(freq_vals(f_s:f_e), mean_values_all(i).pwspectrum_off_normed(ri,f_s:f_e),"Color",[strain_colors(i,:),0.2]);
        hold on;
    end
    loglog(freq_vals(f_s:f_e), mean_spectrum_off_normed(i,f_s:f_e),"Color",strain_colors(i,:),'LineWidth',2);
    hold on;
end
ylim([ymin,ymax]);
ylabel('OFF resp');

%bars of power contribution within each band
xlabelarr = {"< 10Hz", "[10,50] Hz", ">50Hz"};

nexttile;%subplot(3,4,[9:10]);
b = bar(mean_pow_on_normed','w');
hold on;
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_on_normed,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,1);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;

    xvals=b(i).XEndPoints(1);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_on_normed(:,1),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(2);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_on_normed(:,2),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(3);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_on_normed(:,3),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;
end
xticklabels(xlabelarr );
xtickangle(45);
ylabel('ON');


nexttile;%subplot(3,4,[11:12]);
b = bar(mean_pow_off_normed','w');
hold on;
for i=1:ngrps
    nrec = size(mean_values_all(i).pwspectrum_off_normed,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,1);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;

    xvals=b(i).XEndPoints(1);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_off_normed(:,1),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(2);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_off_normed(:,2),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;

    xvals=b(i).XEndPoints(3);
    xvals=repmat(xvals,[1,nrec]);
    xvals2plot=xvals+xdev;
    yvals = reshape(mean_values_all(i).pow_off_normed(:,3),1,[]);
    scatter(xvals2plot,yvals,7,'filled','k',"MarkerFaceColor",strain_colors(i,:));
    hold on;
end
xticklabels(xlabelarr );
xtickangle(45);
ylabel('OFF');
title(t,'normalized power spectrum and relative contribution of frequences');

%% individual traces
for i=1:ngrps
    figure, t=tiledlayout(2,1,"Padding","tight");
    nexttile;
    nrec = size(mean_values_all(i).pwspectrum_on_normed,1);    
    for ri=1:nrec        
        loglog(freq_vals(f_s:f_e), mean_values_all(i).pwspectrum_on_normed(ri,f_s:f_e),"Color",[strain_colors(i,:),0.3]);
        hold on;
    end
    loglog(freq_vals(f_s:f_e), mean_spectrum_on_normed(i,f_s:f_e),"Color",strain_colors(i,:));
    ylabel('ON');
    title(t,legend_name{i},'Interpreter','none');

    nexttile;
    nrec = size(mean_values_all(i).pwspectrum_off_normed,1);    
    for ri=1:nrec        
        loglog(freq_vals(f_s:f_e), mean_values_all(i).pwspectrum_off_normed(ri,f_s:f_e),"Color",[strain_colors(i,:),0.3]);
        hold on;        
    end
    loglog(freq_vals(f_s:f_e),  mean_spectrum_off_normed(i,f_s:f_e),"Color",strain_colors(i,:));
    ylabel('ON');
    title(t,legend_name{i},'Interpreter','none');
end
