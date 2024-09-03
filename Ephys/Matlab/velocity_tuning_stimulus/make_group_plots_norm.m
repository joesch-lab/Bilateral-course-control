% %make population plots
% average of the group responses is normalized, s.t. the max response
% of the group is 1

clear all;
population_folder = 'C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data\Data';

%% define parameters which are required for the data to be selected
sp_freq=0.01;
isHS=1;

%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
% group_labels=[["HS","FlpND"]];
% group_labels=[["HS","FlpND"];["HS","FlpD"]; ["HS","FlpNDxDB331"]]; 
% group_labels=[["VS1_4","FlpD"];["VS1_4","FlpND"];["VS1_4","CantonS"]; ["VS1_4","FlpND DB331"]]; 
% group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","CantonS"]; ["HSE","FlpND DB331"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"];["HSN","CantonS"]; ["HSN","FlpND DB331"]]; 
group_labels=[["H2","FlpND"]];

% option to normalize max response to 1 within each animal before averaging
normalize_to_one = 1;



%% get parameters of the stimulus
if isHS==1
    dir_arr=[0,180];
else
    dir_arr=[90,270];
end


%get all possible speeds for data grouping
[sp_arr,~] = all_speed_spfreq_param_vel_tuning(population_folder);
%get all combinations of speed and directions
[x,y]=meshgrid(sp_arr, dir_arr);
x=x(:); y=y(:);
speed_dir=[x';y']';
nspeedsdir=size(speed_dir,1);


file_ending='_vel_tuning.mat';
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
path_patt=fullfile(population_folder,fullfile('**','vel_tuning_*.mat'));
filelist = dir(path_patt);

%collect data only from the files that match filters and group patterns
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).fly_recordings=[];
    mean_values_all(i).pd_nd=[];    
    mean_values_all(i).files={};
    mean_values_all(i).raw_traces_total=[];
    mean_values_all(i).baseline_traces=[];
    mean_values_all(i).mean_spikes=[];
end

for li=1:length(filelist)    
    load(fullfile(filelist(li).folder,filelist(li).name));
    %exp_data is now available
    if exp_data.isHS~=isHS
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
    speed_dir_i = exp_data.speed_dir;
    nspeedsdir=size(speed_dir,1);
    inds = zeros(size(speed_dir_i,1),1);
    for i=1:length(inds)
        indsi = find(ismember(speed_dir, speed_dir_i(i,:),'rows'));
        inds(i)=indsi;
    end    
    

    %find indices in the sorted speed array 
    mean_i=exp_data.mean_vals;
    mean_traces=NaN(nspeedsdir, 14000);
    baseline_trace=NaN(nspeedsdir, 5000);
    mean_add = NaN(1,nspeedsdir);
    mean_add(inds)=mean_i;
    mean_traces(inds,:)=squeeze(mean(exp_data.raw_traces(:,:,1:14000),2));
    baseline_traces(inds, :)=squeeze(mean(exp_data.baseline_traces(:,:,:),2));

    % spike data
    mean_spikes=exp_data.mean_spikes';    
    mean_add_spikes=mean_spikes;

    pd_nd=mean_add(1:2:end)-mean_add(2:2:end);
    pospdnd=nnz(pd_nd>0);
    if length(pd_nd)-pospdnd>pospdnd
        pd_nd=-pd_nd;
    end
%     if normalize_to_one
%         pd_nd = pd_nd/max(pd_nd);        
% %         mean_add=mean_add/max(mean_add);
%     end

    if isempty(mean_values_all(group_id).fly_recordings)  
        mean_values_all(group_id).fly_recordings=mean_add;
        mean_values_all(group_id).pd_nd=pd_nd;
        mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).raw_traces_total=mean_traces;
        mean_values_all(group_id).baseline_traces=baseline_traces;
        mean_values_all(group_id).mean_spikes=mean_spikes;
    else
        i=size(mean_values_all(group_id).fly_recordings,1)+1;
        mean_values_all(group_id).fly_recordings(i,:)=mean_add;   
        mean_values_all(group_id).pd_nd(i,:)=pd_nd;   
        mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).raw_traces_total=cat(3,mean_values_all(group_id).raw_traces_total, mean_traces);
        mean_values_all(group_id).baseline_traces=cat(3,mean_values_all(group_id).baseline_traces, baseline_traces);
        mean_values_all(group_id).mean_spikes(i,:)=mean_spikes;
    end
end

%make a plot of average values
%find indices of non empty groups
nonemptygrps =[];
for i=1:ngrps
    if ~isempty(mean_values_all(i).fly_recordings)
        nonemptygrps=[nonemptygrps,i];
    end
end
nunneg=length(nonemptygrps);

pdnd_vals=zeros(nunneg,nspeedsdir/2);
stdvals_pdnd=zeros(nunneg,nspeedsdir/2);
meanvals_dir1=zeros(nunneg,nspeedsdir/2);
meanvals_dir2=zeros(nunneg,nspeedsdir/2);
stdvals_dir1=zeros(nunneg,nspeedsdir/2);
stdvals_dir2=zeros(nunneg,nspeedsdir/2);

meanvals_dir1_raw=zeros(nunneg,nspeedsdir/2);
meanvals_dir2_raw=zeros(nunneg,nspeedsdir/2);
stdvals_dir1_raw=zeros(nunneg,nspeedsdir/2);
stdvals_dir2_raw=zeros(nunneg,nspeedsdir/2);

for ii=1:nunneg    
    i=nonemptygrps(ii);
    alldir=mean(mean_values_all(i).fly_recordings,1,'omitNaN');
    meanvals_dir1_raw(ii,:)=alldir(1:2:nspeedsdir);
    meanvals_dir2_raw(ii,:)=alldir(2:2:nspeedsdir);
    alldir=std(mean_values_all(i).fly_recordings,[],1,'omitNaN');
    stdvals_dir1_raw(ii,:)=alldir(1:2:nspeedsdir);
    stdvals_dir2_raw(ii,:)=alldir(2:2:nspeedsdir);

    maxval=max(abs(mean(mean_values_all(i).fly_recordings,1,'omitNaN')));
    nrec=size(mean_values_all(i).fly_recordings,1);
    alldir=mean(mean_values_all(i).fly_recordings,1,'omitNaN')/maxval;
    meanvals_dir1(ii,:)=alldir(1:2:nspeedsdir);
    meanvals_dir2(ii,:)=alldir(2:2:nspeedsdir);
    alldir=std(mean_values_all(i).fly_recordings/maxval,[],1,'omitNaN')/sqrt(nrec);
    stdvals_dir1(ii,:)=alldir(1:2:nspeedsdir);
    stdvals_dir2(ii,:)=alldir(2:2:nspeedsdir);
    maxpdval=max(abs(mean(mean_values_all(i).pd_nd,1,'omitNaN')));
    pdnd_vals(ii,:)=mean(mean_values_all(i).pd_nd,1,'omitNaN')/maxpdval;
    stdvals_pdnd(ii,:)=std(mean_values_all(i).pd_nd/maxpdval,[],1,'omitNaN')/sqrt(nrec);
end

%% save all the data in "figures velocity tuning" folder
all_data = {};
for i=1:size(group_labels, 1)
    name=strcat(group_labels(i,1), '_',group_labels(i,2));
    raw_trace=mean_values_all(i).raw_traces_total;
    s=size(raw_trace);
    all_data.(name).fly_recordings = reshape(permute(raw_trace, [2,1,3]), [], s(3)).';
    baseline_trace=mean_values_all(i).baseline_traces;
    s=size(baseline_trace);
    all_data.(name).baseline_traces = reshape(permute(baseline_trace, [2,1,3]), [], s(3)).';
    all_data.(name).mean_pd = meanvals_dir1_raw(i, :);
    all_data.(name).mean_pd_normed = meanvals_dir1(i, :);
    all_data.(name).mean_nd = meanvals_dir2_raw(i, :);
    all_data.(name).mean_nd_normed = meanvals_dir2(i, :);
    %all_data.name.mean_pdnd = meanvals_dir1(1);
    all_data.(name).mean_pdnd_normed = pdnd_vals(i, :);
    all_data.(name).stderr_pd = stdvals_dir1_raw(i, :);
    all_data.(name).stderr_pd_normed = stdvals_dir1(i, :);
    all_data.(name).stderr_nd = stdvals_dir2_raw(i, :);
    all_data.(name).stderr_nd_normed = stdvals_dir2(i, :);
    %all_data.name.stderr_pdnd = stdvals_dir1(1);
    all_data.(name).stderr_pdnd_normed = stdvals_pdnd(i, :);
    all_data.(name).each_fly_pdnd = mean_values_all(i).pd_nd;
    all_data.(name).each_fly_pd = mean_values_all(i).fly_recordings(:, 1:2:length(mean_values_all(i).fly_recordings));
    all_data.(name).each_fly_nd = mean_values_all(i).fly_recordings(:, 2:2:length(mean_values_all(i).fly_recordings));

    all_data.(name).each_fly_pd_spikes = mean_values_all(i).mean_spikes(:, 1:2:length(mean_values_all(i).mean_spikes));
    all_data.(name).each_fly_nd_spikes = mean_values_all(i).mean_spikes(:, 2:2:length(mean_values_all(i).mean_spikes));
end
filename = fullfile('C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data\Data\FlpD\H2', 'velocity_total_data.mat');
save(filename, 'all_data')
return
%%
legend_name={};
for ii=1:nunneg
    i=nonemptygrps(ii);
    legend_name{ii}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).fly_recordings,1)),')'];
end

strain_colors = brewermap(ngrps,'Set1');
figure, 
subplot(2,1,1);
b = bar(meanvals_dir1');
hold on;
for ii=1:nunneg
    i=nonemptygrps(ii);
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir1(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end

xticklabels(speed_dir(1:2:end,1)');
title('1st direction');
legend(legend_name,'Interpreter','none');

subplot(2,1,2);
b=bar(meanvals_dir2');
hold on;
for ii=1:nunneg
    i=nonemptygrps(ii);
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir2(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end
xticklabels(speed_dir(2:2:end,1)');
title('2nd direction');
% legend(legend_name,'Interpreter','none');

titlestr = ['Velocity tuning for spatial frequency ',num2str(sp_freq)];
sgtitle(titlestr);


figure, 
b = bar(pdnd_vals');
hold on;
for ii=1:nunneg
    i=nonemptygrps(ii);
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_pdnd(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
end

xticklabels(speed_dir(1:2:end,1)');
legend(legend_name,'Interpreter','none');
titlestr = ['PD-ND and velocity tuning for spatial frequency ',num2str(sp_freq)];
title(titlestr);




%% PD and ND bar plots with dots per fly

figure, 
subplot(2,1,1);
b = bar(meanvals_dir1_raw');
ncond = round(size(mean_values_all(nonemptygrps(1)).fly_recordings,2)/2);
hold on;
for ii=1:nunneg
    i=nonemptygrps(ii);
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir1_raw(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
    
    nrec = size(mean_values_all(i).fly_recordings,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,ncond);
    xvals=b(i).XEndPoints;
    xvals=repmat(xvals,[1,nrec]);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;       
    yvals = reshape(mean_values_all(i).fly_recordings(:,1:2:end)',1,[]);
    scatter(xvals,yvals,7,'filled','k');    
    hold on;    
end

xticklabels(speed_dir(1:2:end,1)');
title('1st direction');
legend(legend_name,'Interpreter','none');

subplot(2,1,2);
b=bar(meanvals_dir2_raw');
hold on;
for ii=1:nunneg
    i=nonemptygrps(ii);    
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*stdvals_dir2_raw(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;
             
    nrec = size(mean_values_all(i).fly_recordings,1);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,ncond);
    xvals=b(i).XEndPoints;
    xvals=repmat(xvals,[1,nrec]);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;
    yvals = reshape(mean_values_all(i).fly_recordings(:,2:2:end)',1,[]);
    scatter(xvals,yvals,7,'filled','k');    
    hold on;    
    
end
xticklabels(speed_dir(2:2:end,1)');
title('2nd direction');
% legend(legend_name,'Interpreter','none');

titlestr = ['Velocity tuning for spatial frequency ',num2str(sp_freq)];
sgtitle(titlestr);

%% for all recording within each group plot curves of mean vals to check the quality
for ii=1:nunneg
    i=nonemptygrps(ii);
    figure,
    nrec = size(mean_values_all(i).fly_recordings,1);
    for fi=1:nrec
        subplot(nrec,2,2*(fi-1)+1);
        name=mean_values_all(i).files{fi};
        stind = strfind(name,'vel_tuning_')+length('vel_tuning_');
        enind = strfind(name,'_stim'); enind=enind (end)-1;
        name=name(stind:enind);
        text(0,0,name,'Interpreter','none');
        axis off;
        axis tight;
        subplot(nrec,2,2*fi);
        yvals=mean_values_all(i).fly_recordings(fi,1:2:end);
        plot(yvals,'r'); hold on;
        yvals=mean_values_all(i).fly_recordings(fi,2:2:end);
        plot(yvals,'b'); hold on;
        if fi<nrec
            xticklabels([]);
        else
            xticklabels(speed_dir(2:2:end,1)');;
        end
    end
    sgtitle(legend_name(i))
end









% save('VS14_spf01.mat','mean_values_all');

% 
% 
% % make_pd_nd_raw_plots
% strain_colors = brewermap(ngrps,'Set1');
% alpha=0.3;
% 
% figure,
% bl_len_taken=round(baseline_fraction*mindur_bl);
% 

% %save figure
% filename =['grating_',save_name,file_ending];
% figname_fig=fullfile(population_folder,[filename,'_pdndtrace.fig']);
% figname_png=fullfile(population_folder,[filename,'_pdndtrace.png']);
% savefig(figname_fig);
% saveas(gcf,figname_png,'png');
% 
% 
% 
% % make_average_all_dir_bar_plots();
% %place all averages per direction into one array
% all_av_dir=zeros(ngrps,ndir);
% for gi=1:ngrps
%     all_av_dir(gi,:)=group_vals(gi).av_per_dir;
% end
% 
% 
% xticklabels(dir_vector);
% title_str= ['Directional responses'];
% title(title_str);
% %save figure
% filename =['grating_',save_name,file_ending];
% figname_fig=fullfile(population_folder,[filename,'_av_dir_resp.fig']);
% figname_png=fullfile(population_folder,[filename,'_av_dir_resp.png']);
% savefig(figname_fig);
% saveas(gcf,figname_png,'png');

