% the code for analysing the recordings during ptresentation of grating stimulus
% moving with different velocities. 
% Gratings move either left/right (log files with HS_Vel) or up/down (log files with VS_Vel)
% The scripts computes average response to different speeds

% 04.03.2021
% O.Symonova

parent_folder='\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository - Copy\FlpD\';
parent_folder='\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository - Copy\FlpD\VS\VS1_4\210209\';
% parent_folder='C:\DATA\Data_for_Olga\fly_ephys\FlpD\test';

%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list


%go thru the list, analize each log file with 'scanning_dot' in the title
for li =1:length(filelist)
    if  contains(filelist(li).name,'S_Vel_','IgnoreCase', true) && contains(filelist(li).name,'.log','IgnoreCase', true)  
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        %find the closest pr file
        prfilename = find_pr_file_closest_date(logfile_fullname);
        if prfilename == -1
            warning(['Could not find the matching pr file for ',logfile_fullname]);
            continue;
        end
        folder=filelist(li).folder;
        prfile_fullname=fullfile(folder,prfilename);
        slashinds=strfind(logfile_fullname,filesep);
        dateind=find(slashinds>length(parent_folder),1,'first');
        lastind=min(length(folder),slashinds(dateind));
        filedatestr=folder(length(parent_folder)+1:lastind);
        filedatestr(filedatestr==filesep)=[];
        try
            compute_velocity_responses(parent_folder, prfile_fullname,logfile_fullname, filedatestr);
            close all;
        catch
            warning(['Something went wrong in analysis of ',logfile_fullname]);
        end
    end
end

      
function prfilename = find_pr_file_closest_date(logfilename)
   slashinds=strfind(logfilename,filesep);
   folder=logfilename(1:slashinds(end));
   log_info = dir(logfilename);
   log_datevec=datevec(log_info.date);
   D = dir(folder);
   min_et=3600;
   prfilename=-1;
   for k = 3:length(D) % avoid using the first ones
       currfile = D(k).name; %file name
       if contains(currfile,'.pr','IgnoreCase', true)
           et=abs(etime(datevec(D(k).date),log_datevec));
           if et<min_et
               min_et=et;
               prfilename=currfile;
           end
       end
   end
   if min_et>10
       disp('No matching pr file has been found');
       prfilename=-1;
   end      
end


function  compute_velocity_responses(rootDir, prfullname,logfullname, filedatestr) 
    %% open ephys data
    [Data, Text_header, filenameout, sampling_rate, minV, maxV]=openpr(prfullname,0);
    
    %% list of params:   
   %channel in Data with red frames
    redchannel=2;  
    %value to threshold red frames
    red_threshold=1;
    
    %% init output folder and file names
    %if filedatestr is empty take it as the creation date
    if isempty(filedatestr)
        file_info = dir(prfullname);
        filedatestr = datestr(file_info.date,'yyyymmdd');
    end
    
    disp(['Analysing ',prfullname]);
    [~,prname,~]=fileparts(prfullname);
    
    [filepath,name,ext] = fileparts(logfullname);
    resfolder=fullfile(filepath,'res');
    if ~exist(resfolder)
        mkdir(resfolder);
    end
        
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);
    
    %find start and end frames of moving curtain
    startend=find_start_end_moving_curtain(Data, sampling_rate);
    
    %get param from the log file
    [nrep, tmove, speed_array, speed_index, sp_freq, ifi, isHS] = read_velocity_tuning_stimuli_log(logfullname);
    
    %remove spurious red frames
    tmin=0.25*tmove*sampling_rate;
    indsremove=[];
    for ii=1:size(startend,2)
        if startend(2,ii)-startend(1,ii)<tmin
            indsremove=[indsremove,ii];
        end
    end
    startend(:,indsremove)=[];
    
    %traces between start and end frames
    resp_dur=min(startend(2,:)-startend(1,:));
    
    nspeeds=size(startend,2)/nrep/2; %length(speed_array);
    nruns=nspeeds*2*nrep;
    
    baseline_activity = mean(Data(:,1));
    Data(:,1)=Data(:,1)-baseline_activity;
    ncondintions=nspeeds*2;
    all_traces_info={};
    mean_traces=zeros(nrep,ncondintions,resp_dur);
    minv=inf;
    maxv=-inf;
    for i=1:nruns
        all_traces_info(i).trace=Data(startend(1,i):startend(1,i)+resp_dur-1,1);
        minv=min(minv,min(all_traces_info(i).trace));
        maxv=max(maxv,max(all_traces_info(i).trace));
        ii=ceil(i/ncondintions);
        jj=i-(ii-1)*ncondintions;
        mean_traces(ii,jj,:)=Data(startend(1,i):startend(1,i)+resp_dur-1,1);
        if isHS 
            if mod(i-1,2)==0
                all_traces_info(i).dir=0;
            else
                all_traces_info(i).dir=180;
            end
        else 
            if mod(i,2)==0
                all_traces_info(i).dir=90;
            else
                all_traces_info(i).dir=270;
            end
        end        
        all_traces_info(i).speed = speed_array(speed_index(ceil(i/2)));        
    end
    
    %mean responses
    mean_traces_info={};    
    for i=1:ncondintions
        %get all reps with the same speed and dir values        
        mean_traces_info(i).mean_trace=squeeze(mean(mean_traces(:,i,:),1));
        mean_traces_info(i).mean_val=mean(mean_traces_info(i).mean_trace);
        mean_traces_info(i).speed=all_traces_info(i).speed;
        mean_traces_info(i).dir=all_traces_info(i).dir;
    end
    
    sort_speeds=zeros(ncondintions,3);
    for i=1:ncondintions
        sort_speeds(i,1)=mean_traces_info(i).speed;
        sort_speeds(i,2)=mean_traces_info(i).dir;
        sort_speeds(i,3)=mean_traces_info(i).mean_val;
    end
    sort_speeds=sortrows(sort_speeds,1);  

    speed_labels=sort_speeds(1:2:end,1);
    dir_labels=sort_speeds(1:2,2);
    %   speed_labels_pix_sec=speed_labels/ifi;
    %   speed_labels_deg_sec=round(speed_labels_pix_sec*80/342);
   

    PD=sort_speeds(1:2:end,3);
    ND=sort_speeds(2:2:end,3);
    PD_ND=PD-ND;
    
    figbase = fullfile(resfolder,name); 
    %spatial frequency in degrees
    sp_freq_deg= round((1/sp_freq)*80/342);
    
    %make a plot of PD and ND mean response
    make_figure_PD_and_ND(PD, ND,speed_labels,prname, isHS, sp_freq_deg, figbase);
    
    %make a plot of PD-ND mean response
    make_figure_pd_nd(PD_ND,speed_labels,prname, isHS, sp_freq_deg, figbase);
    
    %make a plot of mean traces per temporal frequency and direction    
    make_figure_mean_trace_per_freq(mean_traces_info, dir_labels, speed_labels, minv, maxv, prname, isHS, sp_freq_deg, figbase);
   
    
    %save variables
    matfile=fullfile(resfolder, [name,'.mat']);
    save(matfile, 'all_traces_info', 'mean_traces_info','sort_speeds', 'speed_labels','baseline_activity','sp_freq_deg','resp_dur','startend','prname','figbase','name', '-nocompression','-v7.3');    
end

function make_figure_PD_and_ND(PD, ND,speed_labels,prname, isHS, sp_freq_deg, figbase)
    %how many different speed conditions
    nspeeds=length(PD);
    
    %plot
    figure, 
    plot(PD,'Marker','o','LineWidth',2); hold on;
    plot(ND,'Marker','o','LineWidth',2); hold on;
    plot(zeros(nspeeds,1),'--');
    xticklabels(speed_labels);
    minV=min([PD', ND']);
    maxV=max([PD', ND']);
    tickvals=[minV,0,maxV];
    yticks(tickvals);
    ystr=num2str(tickvals,'%.3f\n');
    %   ystr=num2str(baseline_activity+tickvals,'%.3f\n');
    yticklabels(ystr);   
    xlabel('Temporal Frequency [Hz]')
    ylabel('Voltage Response');

    %legend
    if isHS
        direction_legend={'rightwards'; 'leftwards';'baseline'};
    else
        direction_legend={'upwards'; 'downwards';'baseline'};
    end
    legend(direction_legend);
    
    %compose the title of the figure
    title_str=prname;
    if isHS
        title_str=[title_str, ', HS velocity'];
    else
        title_str=[title_str, ', VS velocity'];
    end   
    title_str=[title_str,', spatial frequency ',num2str(sp_freq_deg), ' ', char(176),'/cycle'];
    title_str=strrep(title_str,'_','\_');
    title(title_str);   
    
    figsavename=[figbase,'mean_resp.png'];
    saveas(gcf,figsavename,'png');
    figsavename = [figbase,'mean_resp.fig'];
    savefig(figsavename);       
end


function make_figure_pd_nd(PD_ND,speed_labels,prname, isHS, sp_freq_deg, figbase)
    %plot PD-ND
        
    figure,     
    plot(PD_ND,'Marker','o','LineWidth',2); 
    
    xticklabels(speed_labels);
%     minV=min(PD_ND);
%     maxV=max(PD_ND);
%     if minV<0
%         tickvals=[minV,0,maxV];
%     else
%         tickvals=[0,minV,maxV];
%     end
%     yticks(tickvals);
%     ystr=num2str(tickvals,'%.3f\n');
%     %   ystr=num2str(baseline_activity+tickvals,'%.3f\n');
%     yticklabels(ystr);
    
    %   xlabel(['Velocity ',char(176),'/sec']);
    xlabel('Temporal Frequency [Hz]');%['Velocity (cycles/sec)']);
    ylabel('Voltage Response');

    %compose the title of the figure
    title_str=prname;
    if isHS
        title_str=[title_str, ', H velocity, PD-ND'];
    else
        title_str=[title_str, ', V velocity, PD-ND'];
    end
   
    title_str=[title_str,', spatial frequency ',num2str(sp_freq_deg), ' ', char(176),'/cycle'];
    title_str=strrep(title_str,'_','\_');
    title(title_str);   
    
    figsavename=[figbase,'pdnd_mean_resp.png'];
    saveas(gcf,figsavename,'png');
    figsavename = [figbase,'pdnd_mean_resp.fig'];
    savefig(figsavename);      
end


function make_figure_mean_trace_per_freq(mean_traces_info, dir_labels, speed_labels, minv, maxv,  prname, isHS, sp_freq_deg, figbase)
    %how many different speed conditions
    nspeeds = length(unique([mean_traces_info(:).speed]));
    %length of the response
    resp_dur=length(mean_traces_info(1).mean_trace);
    
    %compose the title of the figure
    title_str=prname;
    if isHS
        title_str=[title_str, ', mean traces, H vel'];
    else
        title_str=[title_str, ', mean traces, V vel'];
    end
   
    title_str=[title_str,', spatial frequency ',num2str(sp_freq_deg), ' ', char(176),'/cycle'];
    title_str=strrep(title_str,'_','\_');
      
    
    %assume each subplot is 100x100
    pw=100;
    fw=(nspeeds+2)*pw;
    fh=(2+1.7)*pw;
    sw=1/(nspeeds+2);
    sh=1/3.5;
    figure('Position',[0,0,fw,fh]);
    subplot('Position',[1.75*sw,3*sh,3*sw,0.4*sh]);
    text(0,0,title_str,'FontSize',13,'FontWeight','bold');
    axis off;
    
    for i=1:nspeeds       
        
        key1=find([mean_traces_info(:).dir]==dir_labels(1) & [mean_traces_info(:).speed]==speed_labels(i));
        ax1i=subplot('Position',[sw+(i-1)*sw,1.8*sh,sw,sh]);
        plot(mean_traces_info(key1).mean_trace); 
        ylim([minv,maxv]); xlim([1,resp_dur]); xticks([]);
        if i==1
            yticks([minv,maxv]);
            ylabel('PD')
        else
            yticks([]);
        end
        
        ax2i=subplot('Position',[sw+(i-1)*sw,0.6*sh,sw,sh]);
        key2=find([mean_traces_info(:).dir]==dir_labels(2) & [mean_traces_info(:).speed]==speed_labels(i));
        plot(mean_traces_info(key2).mean_trace);   
        ylim([minv,maxv]); xlim([1,resp_dur]); xticks([]);
        if i==1
            yticks([minv,maxv]);
            ylabel('ND')
        else
            yticks([]);
        end
        xlabel(speed_labels(i));        
    end
    subplot('Position',[2.75*sw,0.2*sh,3*sw,0.4*sh]);
    text(0,0,'Temporal Frequency [Hz]','FontSize',12);
    axis off;
    
    figsavename=[figbase,'_mean_traces.png'];
    saveas(gcf,figsavename,'png');
    figsavename = [figbase,'_mean_traces.fig'];
    savefig(figsavename);   
end

function Data = threshold_red_signal(Data, red_threshold)
    %Data is the Nx2 array, red frames are in the 2nd channel
    rft=find(Data(:,2)>red_threshold); %red frames
    Data(rft,2)=1; %set max values of red to one
    Data(Data(:,2)~=1,2)=0; % set to zero values smaller than 1
end

function startend=find_start_end_moving_curtain(Data, sampling_rate)
   %find the start and end time of moving curtain
   %red frames are shown only during moving part 
   %moving part of the stimulus is separated by >1 sec
   
   %time of red frames
   rft=find(Data(:,2));     
   %find the time difference between two consecutive red frames
   drf=rft(2:end)-rft(1:end-1);
   gap_min= 0.5*sampling_rate;
   
   %find big gaps in data (gap of 0.5sec is hardcoded in stimulus) 
   gap_rf_ind = find(drf>gap_min);
   lastind = [rft(gap_rf_ind); rft(end)];
   firstind= [rft(1); rft(gap_rf_ind+1)];
   
   startend=[firstind';lastind'];
   
%    figure, plot(Data);
%    plot_pair_frames(firstind,lastind);   
end


function plot_pair_frames(timest,timeen)
    %assume the figure with the ephys activity is the current
    hold on 
    plot([timest, timest],[0,1], 'm'); %start frame
    hold on
    plot([timeen, timeen],[0,1], 'k'); %end
    hold off
end




 


function [nrep, tmove, speed_array, speed_indices, sp_freq, ifi, isHS] = read_velocity_tuning_stimuli_log(logfile)

    % read the file with stimuli parameters
    f=fopen(logfile);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'_Vel'))
        disp('The log file does not contain information about velocity selectivity stimuli.');
        fclose(f);
        return;
    end
    
    if isempty(strfind(tline,'HS'))
        isHS=0;
    else
        isHS=1;
    end
    
    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    tend=datevec(tstr,formatIn);
    
    tline=fgetl(f);
    t_str=tline(length('Time per grating: ')+1:end);     
    tmove=uint32(str2double(t_str)); 
    
    tline=fgetl(f);
    ifistr=tline(length('Flip rate: ')+1:end);    
    ifi=str2double(ifistr);
    
    tline=fgetl(f);
    nrepstr=tline(length('Repetitions: ')+1:end);    
    nrep=str2num(nrepstr); 

    speed_array=[];    
    tline=fgetl(f);
    while ~isempty(strfind(tline,'Speed'))
        tstr=tline(length('Speed: ')+1:end);    
        speedi=str2double(tstr); 
        speed_array=[speed_array,speedi];
         tline=fgetl(f);
    end
        
    %Baseline recording: skipped
    
    tline=fgetl(f);
    sf_str=tline(length('Spatial frequency: ')+1:end);     
    sp_freq=str2double(sf_str); 
    
    speed_indices=[];
    tline=fgetl(f);
    while ~isempty(strfind(tline,'Random_Number'))
        tstr=tline(length('Random_Number: ')+1:end);    
        speedi=str2num(tstr); 
        speed_indices=[speed_indices,speedi];
         tline=fgetl(f);
    end
    
    fclose(f);
end