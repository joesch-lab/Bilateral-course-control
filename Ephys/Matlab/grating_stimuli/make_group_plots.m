% %make population plots

% average of the group responses is normalized, s.t. the max response
% of the group is 1

clear all;
population_folder = 'C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data\Data';


%% define parameters which are required for the data to be selected
sp_freq=0.01;
dx=5;

%can define any group to plot, the labels of the groups 
% are defined within [] and separated by ;
% group_labels=[["HS","FlpND"];["HS","FlpD"]; ["HS", "FlpNDxDB331"]; ["HS", "CantonS"]]; 
% group_labels=[["VS1_4","FlpD"];["VS1_4","FlpND"];["VS1_4","CantonS"]; ["VS1_4","FlpND DB331"]]; 
% group_labels=[["HSE","FlpD"];["HSE","FlpND"];["HSE","CantonS"]; ["HSE","FlpND DB331"]]; 
% group_labels=[["HSN","FlpD"];["HSN","FlpND"];["HSN","CantonS"]; ["HSN","FlpND DB331"]]; 
% group_labels=[["H2","FlpD"];["H2","FlpND"]];
group_labels=[["H2","FlpND"]];

ngrps= size(group_labels,1);
for ci=1:ngrps
    file_patt_i='';
    for cii=1:length(group_labels(ci,:))-1
        file_patt_i=strcat(file_patt_i,strcat(group_labels(ci,cii),'_'));        
    end
    file_patt_i = strcat(file_patt_i,group_labels(ci,end));
    file_patt{ci}=char(file_patt_i);
end


% option to normalize max response to 1 within each animal before averaging
normalize_to_one = 1;
%determine H or V cells to analyze from group labels
if ~isempty(strfind(group_labels(1),'H'))
    isHS=1;
else
    isHS=0;
end

%% get preferred and null directions
if isHS==1
    dir_arr=[0,180];
else
    dir_arr=[90,270];
end


file_pattern='grating_folder_*.mat';


%load all files that start with 'grating_folder_' and contain labels
%specified in one of the groups
%get list of all files and folders 
path_patt=fullfile(population_folder,fullfile('**',file_pattern));
filelist = dir(path_patt);

%collect data only from the files that match filters and group patterns
mean_values_all={};
for i=1:ngrps
    mean_values_all(i).name=group_labels(i,:);
    mean_values_all(i).fly_recordings=[];
    mean_values_all(i).bl=[];
    mean_values_all(i).dir_responses=[];
    mean_values_all(i).mean_spikes=[];
    mean_values_all(i).files={};
end

%crop all recordins to this lenght:
sampling_rate=10000;
minResp=0.5*sampling_rate;
minbl=0.5*sampling_rate;

dir_vector=0:45:359;
ndir=length(dir_vector);

for li=1:length(filelist)    
   %find the matching group, otherwise skip    
    filepath=filelist(li).folder;
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
    %if stim params do not match, skip the recording
    if exp_data.stimparam.speed~=dx
        continue;
    end
    if exp_data.stimparam.sp_freq~=sp_freq
        continue;
    end

    
    %find indices in the sorted speed array 
    mean_i=mean(exp_data.traces_av(:,1:minResp),2);
%     if normalize_to_one
%         mean_i= mean_i/max(abs(mean_i));        
%     end
%     
    if isempty(mean_values_all(group_id).fly_recordings)          
        mean_values_all(group_id).fly_recordings=exp_data.traces_av(:,1:minResp);
        mean_values_all(group_id).bl=exp_data.baseline_av(:,end-minbl+1:end);
        mean_values_all(group_id).dir_responses=mean_i;    
        mean_values_all(group_id).files={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).mean_spikes=exp_data.meanspike_rate;
    else        
        mean_values_all(group_id).fly_recordings=cat(3,mean_values_all(group_id).fly_recordings,exp_data.traces_av(:,1:minResp));
        mean_values_all(group_id).bl=cat(3,mean_values_all(group_id).bl,exp_data.baseline_av(:,end-minbl+1:end));
        mean_values_all(group_id).dir_responses=cat(2,mean_values_all(group_id).dir_responses,mean_i);        
        i=size(mean_values_all(group_id).fly_recordings,3);
        mean_values_all(group_id).files(i)={fullfile(filelist(li).folder,filelist(li).name)};
        mean_values_all(group_id).mean_spikes(:, i)=exp_data.meanspike_rate;
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
alldur=minbl+minResp;

dir_vals=zeros(nunneg,ndir);
dir_vals_std=zeros(nunneg,ndir);
dir_vals_raw=zeros(nunneg,ndir);
dir_vals_std_raw=zeros(nunneg,ndir);
meantraces=zeros(nunneg,ndir,alldur);
stdtraces=zeros(nunneg,ndir,alldur);

for ii=1:nunneg
    i=nonemptygrps(ii);
    nrec = size(mean_values_all(i).fly_recordings,3);

    allresp=mean(mean_values_all(i).fly_recordings,3,'omitNaN');
    allbl=mean(mean_values_all(i).bl,3,'omitNaN');
    allresp_std=std(mean_values_all(i).fly_recordings,[],3,'omitNaN');
    allbl_std=std(mean_values_all(i).bl,[],3,'omitNaN');

    meantraces(ii,:,:)=cat(2,allbl,allresp);
    stdtraces(ii,:,:)=cat(2,allbl_std,allresp_std);
    max_ii=max(abs(mean(mean_values_all(i).dir_responses,2)));
    dir_vals(ii,:)=mean(mean_values_all(i).dir_responses,2)/max_ii;  
    dir_vals_std(ii,:)=std(mean_values_all(i).dir_responses/max_ii,[],2)/sqrt(nrec);
    dir_vals_raw(ii,:)=mean(mean_values_all(i).dir_responses,2); 
    dir_vals_std_raw(ii,:)=std(mean_values_all(i).dir_responses,[],2)/sqrt(nrec); 
end
legend_name={};
for ii=1:nunneg
    i=nonemptygrps(ii);
    legend_name{ii}=[char(strjoin(group_labels(i,:),' ')),' (',num2str(size(mean_values_all(i).fly_recordings,3)),')'];
end

strain_colors = brewermap(ngrps,'Set1');
%% save all the data to folder "figures gratings with checkers"
all_data = {};
for i=1:size(group_labels, 1)
    name=strcat(group_labels(i,1), '_',group_labels(i,2));
    raw_trace=mean_values_all(i).fly_recordings;
    s=size(raw_trace);
    all_data.(name).fly_recordings = reshape(permute(raw_trace, [2,1,3]), [], s(3)).';
    baseline_trace=mean_values_all(i).bl;
    s=size(baseline_trace);
    all_data.(name).baseline_traces = reshape(permute(baseline_trace, [2,1,3]), [], s(3)).';
    all_data.(name).each_fly_resp = squeeze(mean_values_all(i).dir_responses);

    all_data.(name).mean_spikes = squeeze(mean_values_all(i).mean_spikes);
end
filename = fullfile('C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data\Data', 'gratings_total_data.mat');
save(filename, 'all_data')
return 
%% bar plot with error bars
figure, 
b = bar(dir_vals','grouped');
xticklabels(dir_vector);

hold on;
for ii=1:nunneg
    i=nonemptygrps(ii);    
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*dir_vals_std(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on;  
end
legend(legend_name,'Interpreter','none');
hold on; %add dots per recording
titlestr = ['Directional responses for speed ',num2str(dx),' and spatial frequency ',num2str(sp_freq)];
title(titlestr);

%% bar plot of raw responses with dots per fly
figure, 
b = bar(dir_vals_raw','grouped');
xticklabels(dir_vector);

hold on;
for ii=1:nunneg
    i=nonemptygrps(ii); 
    %error bar
    xvals=b(i).XEndPoints;
    xvals=[xvals;xvals];    
    yvals=b(i).YEndPoints;
    stderr_signed=sign(yvals).*dir_vals_std_raw(i,:);    
    yvalsup=yvals+ stderr_signed;
    yvals=[yvals; yvalsup];
    plot(xvals,yvals,'-k_');
    b(i).FaceColor=strain_colors(i,:);
    hold on; 
    %dots
    nrec = size(mean_values_all(i).fly_recordings,3);
    xdev =rand(1,nrec);
    xdev =repelem(xdev,ndir);
    xvals=b(i).XEndPoints;
    xvals=repmat(xvals,[1,nrec]);
    barw=b(i).BarWidth*0.2;
    xdev=xdev*barw-barw/2;
    xvals=xvals+xdev;    
    yvals = reshape(mean_values_all(i).dir_responses,1,[]);
    scatter(xvals,yvals,7,'filled','k');    
    hold on;    
end
legend(legend_name,'Interpreter','none');
hold on; %add dots per recording
titlestr = ['Directional responses for speed ',num2str(dx),' and spatial frequency ',num2str(sp_freq)];
title(titlestr);

%% mean traces of the group 
figure,
for dii=1:ndir+1    
    subplot(ndir+1,1,dii);
    if dii==1
        axis off;
        strtex='Traces of ';
        for ii=1:nunneg
            i=nonemptygrps(ii);
            grname=strrep(file_patt{i},'_',' ');
            strtex=[strtex, ' \color[rgb]{' sprintf('%1.2f,%1.2f,%1.2f', strain_colors(ii,:)) '} ' grname];
        end
        text(0.3,0.3,strtex,'Interpreter','tex');
    else
        di=dii-1;
        for ii=1:nunneg
            i=nonemptygrps(ii);
            meantracei=squeeze(meantraces(i,di,:))';
            stdtracei=squeeze(stdtraces(i,di,:))';
            stdshade_mean_std(meantracei,stdtracei,0.2,strain_colors(ii,:));
            ylabel(dir_vector(di));
            hold on;
        end
    end
end


%% for each group make a plot of individual traces for debugging
for ii=1:nunneg
    i=nonemptygrps(ii);
    figure;
    t = tiledlayout(ndir+1,1,'TileSpacing','Compact');
    legname={};
    nrec=size(mean_values_all(i).fly_recordings,3);

    for dii=1:ndir
        nexttile
        for ri=1:nrec
            if dii==1
                filerec=mean_values_all(i).files{ri};
                subscr=strfind(filerec,'_');
                filerec=filerec(subscr(end-4)+1:subscr(end-1)-1);
                legname{ri}=filerec;
            end
            tri=[squeeze(mean_values_all(i).bl(dii,:,ri)),squeeze(mean_values_all(i).fly_recordings(dii,:,ri))];
            plot((1:numel(tri))/sampling_rate,tri);
            hold on;
        end 
        ylabel(dir_vector(dii))
    end    
    hl = legend(legname,"Interpreter","none","Location","east","NumColumns",3);
    set(hl,'position',[0.15 0.05 .7 .1])
end

    
