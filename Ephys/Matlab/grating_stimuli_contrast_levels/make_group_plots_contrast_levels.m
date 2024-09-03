% %make population plots
% responses of the groups are normalized with max to one;

clear all;
population_folder = 'C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data\Data';


%% define parameters which are required for the data to be selected
sp_freq=0.01;
dx=1;

%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
% group_labels=[["VS","FlpD"];["VS","FlpND"];["VS","CantonS"]; ["VS","FlpND DB331"]]; 
% group_labels=[["HS","FlpND"];["HS","FlpD"]; ["HS","FlpNDxDB331"]];
% group_labels=[["HS","FlpD"];["HS","FlpND"];["HS","CantonS"]; ["HS","FlpND DB331"]]; 
% group_labels=[["VS1_4","FlpD"];["VS1_4","FlpND"];["VS1_4","CantonS"]; ["VS1_4","FlpND DB331"]]; 
% group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","CantonS"]; ["HSE","FlpND DB331"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"];["HSN","CantonS"]; ["HSN","FlpND DB331"]]; 
% group_labels=[["H2","FlpD"]; ["H2","FlpND"]];
group_labels=[["H2","FlpND"]];



%% get parameters of the stimulus
file_patt='grating_contrast_levels_*.mat';
ngrps= size(group_labels,1);


%load all file that start with 'grating_contrast_levels_' and contain labels
%specified in the group
%get list of all files and folders 
path_patt=fullfile(population_folder,fullfile('**','grating_contrast_levels_*.mat'));
filelist = dir(path_patt);

%collect data only from the files that match filters and group patterns

mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).mean_contrast=[];
    mean_values_all(i).pd_nd=[]; 
    mean_values_all(i).files={};
    mean_values_all(i).raw_traces=[];
    mean_values_all(i).raw_traces_total=[];
    mean_values_all(i).baseline_traces=[];
    mean_values_all(i).mean_spikes=[];
%     mean_values_all(i).raw_spikes_total=[];
%     mean_values_all(i).baseline_spikes=[];
end

contrast_dir=[];
nspeedsdir=0;
sampling_rate=10000;
minrawtrace=0.5*sampling_rate;

for li=1:length(filelist)    
    load(fullfile(filelist(li).folder,filelist(li).name));
    %exp_data is now available
    if exp_data.speed~=dx
        continue;
    end
    if exp_data.sp_freq~=sp_freq
        continue;
    end
    
    filepath=exp_data.folder;
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
    %add mean and std values to the group
    contrast_dir_i = exp_data.allconfigs;

    
    %all contrast levels sequense should be identical among recordings
    if nspeedsdir==0
        contrast_dir=contrast_dir_i;
        nspeedsdir=size(contrast_dir,1);
        nspeeds_1dir=nspeedsdir/2;
    end        

    %find indices in the sorted speed array
    mean_i=exp_data.mean_vals';    
    mean_add=mean_i;
    %spike data
    mean_spikes=exp_data.mean_spikes';    
    mean_add_spikes=mean_spikes;
    
    pd_nd=mean_add(1:nspeeds_1dir)-mean_add(nspeeds_1dir+1:end);
    pospdnd=nnz(pd_nd>0);
    if length(pd_nd)-pospdnd>pospdnd
        pd_nd=-pd_nd;
    end
        
    mean_raw_traces=squeeze(mean(exp_data.raw_traces(:,:,1:minrawtrace),2));
    raw_traces_total=squeeze(mean(exp_data.raw_traces(:,:,:),2));
    baseline_traces=squeeze(mean(exp_data.baseline_traces(:,:,:),2));
    if isempty(mean_values_all(group_id).mean_contrast)  
        mean_values_all(group_id).mean_contrast=mean_add;
        mean_values_all(group_id).pd_nd=pd_nd;
        mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).raw_traces=mean_raw_traces;
        mean_values_all(group_id).raw_traces_total=raw_traces_total;
        mean_values_all(group_id).baseline_traces=baseline_traces;
        mean_values_all(group_id).mean_spikes=mean_spikes;
    else
        i=size(mean_values_all(group_id).mean_contrast,1)+1;
        mean_values_all(group_id).mean_contrast(i,:)=mean_add;   
        mean_values_all(group_id).pd_nd(i,:)=pd_nd;   
        mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).raw_traces=cat(3,mean_values_all(group_id).raw_traces,mean_raw_traces);
        mean_values_all(group_id).raw_traces_total=cat(3,mean_values_all(group_id).raw_traces_total, raw_traces_total);
        mean_values_all(group_id).baseline_traces=cat(3,mean_values_all(group_id).baseline_traces, baseline_traces);
        mean_values_all(group_id).mean_spikes(i,:)=mean_spikes;
    end
end

%find indices of non empty groups
nonemptygrps =[];
for i=1:ngrps
    if ~isempty(mean_values_all(i).mean_contrast)
        nonemptygrps=[nonemptygrps,i];
    end
end
nunneg=length(nonemptygrps);

%% compute averages accross multiple cells per each group
meanvals_dir1=zeros(nunneg,nspeedsdir/2);
meanvals_dir2=zeros(nunneg,nspeedsdir/2);
stdvals_dir1=zeros(nunneg,nspeedsdir/2);
stdvals_dir2=zeros(nunneg,nspeedsdir/2);

pdnd_vals_normed=zeros(nunneg,nspeedsdir/2);
stdvals_pdnd_normed=zeros(nunneg,nspeedsdir/2);
meanvals_dir1_normed=zeros(nunneg,nspeedsdir/2);
meanvals_dir2_normed=zeros(nunneg,nspeedsdir/2);
stdvals_dir1_normed=zeros(nunneg,nspeedsdir/2);
stdvals_dir2_normed=zeros(nunneg,nspeedsdir/2);

for ii=1:nunneg
    i=nonemptygrps(ii);    
    alldir=mean(mean_values_all(i).mean_contrast,1,'omitNaN');
    maxval=max(abs(alldir));

    stdvals = std(mean_values_all(i).mean_contrast,[],1,"omitnan");
    meanvals_dir1(ii,:)=alldir(1:nspeedsdir/2);
    meanvals_dir2(ii,:)=alldir(nspeedsdir/2+1:nspeedsdir);
    stdvals_dir1(ii,:)=stdvals(1:nspeedsdir/2);
    stdvals_dir2(ii,:)=stdvals(nspeedsdir/2+1:nspeedsdir);
    

    meanvals_dir1_normed(ii,:)=alldir(1:nspeedsdir/2)/maxval;
    meanvals_dir2_normed(ii,:)=alldir(nspeedsdir/2+1:nspeedsdir)/maxval;
    stdvals=std(mean_values_all(i).mean_contrast/maxval,[],1,'omitNaN');
    nrec=size(mean_values_all(i).mean_contrast,1);
    stdvals_dir1_normed(ii,:)=stdvals(1:nspeedsdir/2)/sqrt(nrec);
    stdvals_dir2_normed(ii,:)=stdvals(nspeedsdir/2+1:nspeedsdir)/sqrt(nrec);
    
    maxpdval=max(abs(mean(mean_values_all(i).pd_nd,1,'omitNaN')));
    pdnd_vals_normed(ii,:)=mean(mean_values_all(i).pd_nd,1,'omitNaN')/maxpdval;
    stdvals_pdnd_normed(ii,:)=std(mean_values_all(i).pd_nd/maxpdval,[],1,'omitNaN')/sqrt(nrec);
end
%% save all the data in "figures contrast level" folder
all_data = {};
for i=1:size(group_labels, 1)
    name=strcat(group_labels(i,1), '_',group_labels(i,2));
    raw_trace=mean_values_all(i).raw_traces_total;
    s=size(raw_trace);
    all_data.(name).fly_recordings = reshape(permute(raw_trace, [2,1,3]), [], s(3)).';
    baseline_trace=mean_values_all(i).baseline_traces;
    s=size(baseline_trace);
    all_data.(name).baseline_traces = reshape(permute(baseline_trace, [2,1,3]), [], s(3)).';
    all_data.(name).mean_pd = meanvals_dir1(i, :);
    all_data.(name).mean_pd_normed = meanvals_dir1_normed(i, :);
    all_data.(name).mean_nd = meanvals_dir2(i, :);
    all_data.(name).mean_nd_normed = meanvals_dir2_normed(i, :);
    %all_data.name.mean_pdnd = meanvals_dir1(1);
    all_data.(name).mean_pdnd_normed = pdnd_vals_normed(i, :);
    all_data.(name).stderr_pd = stdvals_dir1(i, :);
    all_data.(name).stderr_pd_normed = stdvals_dir1_normed(i, :);
    all_data.(name).stderr_nd = stdvals_dir2(i, :);
    all_data.(name).stderr_nd_normed = stdvals_dir2_normed(i, :);
    %all_data.name.stderr_pdnd = stdvals_dir1(1);
    all_data.(name).stderr_pdnd_normed = stdvals_pdnd_normed(i, :);
    all_data.(name).each_fly_pdnd = mean_values_all(i).pd_nd;
    all_data.(name).each_fly_pd = mean_values_all(i).mean_contrast(:, 1:5);
    all_data.(name).each_fly_nd = mean_values_all(i).mean_contrast(:, 6:10);

    all_data.(name).each_fly_pd_spikes = mean_values_all(i).mean_spikes(:, 1:5);
    all_data.(name).each_fly_nd_spikes = mean_values_all(i).mean_spikes(:, 6:10);
end
filename = fullfile('C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data\Data', 'contrast_total_data.mat');
save(filename, 'all_data')
return
%% plot normalized responses per contrast level for each direction
legend_name={};
for ii=1:nunneg
    i=nonemptygrps(ii);
    legend_name{ii}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).mean_contrast,1)),')'];
end

strain_colors = brewermap(ngrps,'Set1');
figure, 
subplot(2,1,1);
b = bar(meanvals_dir1_normed');
hold on;
for i=1:nunneg    
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir1_normed(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end

xticklabels(["high"',"", "", "", "low"]);
title('1st direction');
legend(legend_name,'Interpreter','none');

subplot(2,1,2);
b=bar(meanvals_dir2_normed');
hold on;
for i=1:nunneg    
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir2_normed(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end
xticklabels(["high"',"", "", "", "low"]);
title('2nd direction');
% legend(legend_name,'Interpreter','none');

titlestr = ['Response to different contrasts (sp freq: ',num2str(sp_freq), ', dx: ', num2str(dx),')'];
sgtitle(titlestr);

%% normalized pd-nd responses
figure, 
b = bar(pdnd_vals_normed');
hold on;
for i=1:nunneg    
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_pdnd_normed(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end

xticklabels(["high"',"", "", "", "low"]);
legend(legend_name,'Interpreter','none');
titlestr = ['PD-ND and different contrasts (sp freq: ',num2str(sp_freq), ', dx: ', num2str(dx),')'];
title(titlestr);

%% bar plot with dots not normalized for data investigation
figure, 
subplot(2,1,1);
b = bar(meanvals_dir1');
hold on;
%add dots per recording
for ii=1:nunneg
    i=nonemptygrps(ii);   
    %error bar
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir1(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
    %dots
    nrec = size(mean_values_all(i).mean_contrast,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,nspeedsdir/2);
    xvals=b(i).XEndPoints;
    xvals=repmat(xvals,[1,nrec]);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;    
    yvals = reshape(mean_values_all(i).mean_contrast(:,1:nspeedsdir/2)',1,[]);
    scatter(xvals,yvals,7,'filled','k');    
    hold on;    
end

xticklabels(["high"',"", "", "", "low"]);
title('1st direction');
legend(legend_name,'Interpreter','none');

subplot(2,1,2);
b=bar(meanvals_dir2');
hold on;
%add dots per recording
for ii=1:nunneg
    i=nonemptygrps(ii);   
    %error bar
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir2(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
    %dots
    nrec = size(mean_values_all(i).mean_contrast,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,nspeedsdir/2);
    xvals=b(i).XEndPoints;
    xvals=repmat(xvals,[1,nrec]);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;    
    yvals = reshape(mean_values_all(i).mean_contrast(:,nspeedsdir/2+1:end)',1,[]);
    scatter(xvals,yvals,7,'filled','k');    
    hold on;    
end
xticklabels(["high"',"", "", "", "low"]);
title('2nd direction');
% legend(legend_name,'Interpreter','none');

titlestr = ['Response to different contrasts (sp freq: ',num2str(sp_freq), ', dx: ', num2str(dx),')'];
sgtitle(titlestr);


%% for each group make a plot of individual traces for debugging
for ii=1:nunneg
    i=nonemptygrps(ii);
    figure;
    t = tiledlayout(nspeedsdir+1,1,'TileSpacing','Compact');
    legname={};
    nrec=size(mean_values_all(i).raw_traces,3);

    for dii=1:nspeedsdir
        nexttile
        for ri=1:nrec
            if dii==1
                filerec=mean_values_all(i).files{ri};
                subscr=strfind(filerec,'_');
                filerec=filerec(subscr(end-4)+1:subscr(end-1)-1);
                legname{ri}=filerec;
            end
            tri=squeeze(mean_values_all(i).raw_traces(dii,:,ri));
            plot((1:numel(tri))/sampling_rate,tri);
            hold on;
        end 
        ylabel([num2str(contrast_dir(dii,1)),' ',num2str(contrast_dir(dii,2))]);
    end    
    hl = legend(legname,"Interpreter","none","Location","east","NumColumns",3);
    set(hl,'position',[0.15 0.05 .7 .1])
end
