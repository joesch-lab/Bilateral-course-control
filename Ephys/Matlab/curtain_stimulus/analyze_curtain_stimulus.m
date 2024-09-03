% the code for analysing the recordings during curtain stimulus
% the stimulus is similar to the grating stimulus. 
% Instead of the the grating moving, either white or
% black part of the curtain expands and hence turning the whole screen into
% white or black.

% 02.03.2021
% O.Symonova

% parent_folder='\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository - Copy\';
parent_folder='C:\Users\rsatapat\Documents\Victoria\Recordings for the paper\Data';
% parent_folder='\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository - Copy\FlpD\HS\HSS\200205\curtain_stimulus\';
% C:\DATA\Data_for_Olga\fly_ephys\FlpND\201202\copy\curtain_stimulus

%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list


%go thru the list, analize each log file with 'scanning_dot' in the title
for li =1:length(filelist)
    if  contains(filelist(li).name,'curtain','IgnoreCase', true) && contains(filelist(li).name,'.log','IgnoreCase', true)  
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
            compute_curtain_responses(parent_folder, prfile_fullname,logfile_fullname, filedatestr);
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
   if min_et>180
       disp('No matching pr file has been found');
       prfilename=-1;
   end      
end


function  compute_curtain_responses(rootDir, prfullname,logfullname, filedatestr) 
    %% open ephys data
    [Data, Text_header, filenameout, sampling_rate, minV, maxV]=openpr(prfullname,0);
    
    %% list of params:   
   %channel in Data with red frames
    redchannel=2;  
    %value to threshold red frames
    red_threshold=4;
    
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
    [nrep, nframes, speed, curtain_bw, sp_freq, ifi] = read_curtain_stimuli_log(logfullname);
    
    %remove spurious red frames
    tmin=0.25*double(nframes)*ifi*sampling_rate;
    indsremove=[];
    for ii=1:size(startend,2)
        if startend(2,ii)-startend(1,ii)<tmin
            indsremove=[indsremove,ii];
        end
    end
    startend(:,indsremove)=[];
    
    
    %traces between start and end frames 
    real_resp_dur=round(nframes*sampling_rate*ifi);
    resp_dur=real_resp_dur+4000; %2000 frame before and 2000 frames after
    nruns=nrep*8;
    
    baseline_activity = mean(Data(:,1));
    Data(:,1)=Data(:,1)-baseline_activity;
    
    all_traces=zeros(nruns, resp_dur);
    for i=1:nruns
        all_traces(i,:)=Data(startend(1,i)-2000:startend(1,i)-2000+resp_dur-1,1);
    end
    
    %mean responses
    mean_traces=zeros(8, resp_dur);
    for i=1:8
        alltraces_i=all_traces(i:8:nruns,:);
        mean_traces(i,:)=mean(alltraces_i);
    end
    
    %average of the response in the whole period
    % remove first 2000 and last 2000 frames before computing the total
    % average
    av_total=mean(mean_traces(:, 2000:2000+real_resp_dur),2);     
    %from stimulus file, then flipped as curtaing closes on the screen
    direction_vector = [180,135,90,45,0,315,270,225];%0:45:359;
    
    %compose the title of the figure
    speed_pix_sec=speed/ifi;
    %the stimulus spans around 80 degrees of visual field
    speed_deg=round(speed_pix_sec*80/342);
    %spatial frequency in degrees
    sp_freq_deg= round((1/sp_freq)*80/342);
    
    title_str=[prname, ', spatial frequency: ',num2str(sp_freq_deg), ', speed: ',num2str(speed_deg)];
    if curtain_bw ==1 %white curtain %on response
        title_str=[title_str,', ON response'];
    else
        title_str=[title_str,', OFF response'];
    end
    title_str=strrep(title_str,'_','\_');
    
    matfile=fullfile(resfolder, [name,'.mat']);
    figbase = fullfile(resfolder,name); 
        
    %make polar plot of direction responses
    make_polar_plot_direction_selectivity(av_total, mean_traces,baseline_activity, direction_vector, title_str, figbase);
    %save variables
    save(matfile, 'direction_vector', 'av_total','mean_traces', 'all_traces','baseline_activity','sp_freq_deg','speed_deg', 'curtain_bw','real_resp_dur','startend','prname','figbase','name', '-nocompression','-v7.3');    
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

function make_polar_plot_direction_selectivity(av_total, mean_traces, baseline_activity, direction_vector, title_str, figbase)
    figure('Position',[0,0,700,2000]);
    
    %title
    subplot('Position',[0.2,0.95,0.6,0.05]);
    text(0,0,title_str, 'fontsize', 13, 'FontWeight','bold');
    axis off;

    %polar plot
    subplot('Position',[0.3,0.5,0.4,0.4]);
    alphas=direction_vector*pi/180;
    alphas_plot=[alphas,alphas(1)];
    %translate all the plots, so that the min value is positive
    minval=min(av_total);
    if minval<0
        av_total_plot=av_total+abs(minval);
        baseline_plot=abs(minval);
    else
        av_total_plot=av_total;
        baseline_plot=0;
    end
    av_total_plot=[av_total_plot',av_total_plot(1)];
    polarplot(alphas_plot,av_total_plot, 'LineWidth',2);
    hold on;
    polarplot(alphas_plot,ones(9,1)*baseline_plot,'--', 'LineWidth',1);
    thetaticks(sort(direction_vector));
    rticks([max(av_total_plot)]);
%     rticks([abs(baseline_plot),max(av_total_plot)]);
%     rticklabels([baseline_activity,baseline_activity+max(av_total)]);%     
    rticklabels([max(av_total)]);
    
    %plot all mean traces together
    subplot('Position',[0.1,0.05,0.8,0.4]);
    for ni=1:8       
        plot(mean_traces(ni,:),'LineWidth',2); hold on;
    end 
    nvals=size(mean_traces,2);
    minval=min(mean_traces(:));
    maxval=max(mean_traces(:));
    yticks([minval,0,maxval]);    
%     ystr=num2str([minval+baseline_activity,baseline_activity,maxval+baseline_activity],'%.3f\n');
    ystr=num2str([minval,0,maxval],'%.3f\n');
    yticklabels(ystr);
    xticks(0:1000:nvals);
    xticklabels((0:1000:nvals)/10000);
    legend(num2str(direction_vector(:)));
    ylabel('Voltage Response'); 

    figsavename=[figbase,'polar_DS.png'];
    saveas(gcf,figsavename,'png');
    figsavename = [figbase,'polar_DS.fig'];
    savefig(figsavename);
end


function make_figure_linked_rep_activity(activity_rep,resfolder)
    %% plot repetitions in separate subplots
    nrep=size(activity_rep,1);
    sphandles=[];
    min_a = min(activity_rep(:));
    max_a = max(activity_rep(:));
    figure,
    for ni=1:nrep       
        sphandi = subplot(nrep,1,ni);
        sphandles=[sphandles,sphandi];
        plot(activity_rep(ni,:));
        ylim([min_a,max_a]);
        xticks([]); yticks([]);       
    end   
   linkaxes(sphandles,'x');
   title(sphandles(1),'Activity during each repetition');
   figsavename=fullfile(resfolder,'rep_activity.png');
   saveas(gcf,figsavename,'png');
   figsavename = fullfile(resfolder,'rep_activity.fig');
   savefig(figsavename);
   %%
end
 
function make_figure_rep_activity_together(activity_rep,resfolder)
   %% plot all reps on one plot and then add average in bold   
   nrep=size(activity_rep,1);
   figure,
   for ni=1:nrep
       plot(activity_rep(ni,:)); hold on;
   end
   av_activity=mean(activity_rep,1);
   plot(av_activity,'k','LineWidth',2);
   title('Activity during all repetitions');
   figsavename=fullfile(resfolder,'all_reps_mean_activity.png');
   saveas(gcf,figsavename,'png');
   figsavename = fullfile(resfolder,'all_reps_mean_activity.fig');
   savefig(figsavename);
   %%
end

function [nrep, nframes,  speed, curtain_bw, sp_freq, ifi] = read_curtain_stimuli_log(logfile)

    % read the file with stimuli parameters
    f=fopen(logfile);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'curtain_stimuli'))
        disp('The log file does not contain information about curtain stimuli.');
        fclose(f);
        return;
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
    fr_str=tline(length('Frames needed to close the curtain: ')+1:end);     
    nframes=uint32(str2num(fr_str)); 
    
    tline=fgetl(f);
    ifistr=tline(length('Flip rate: ')+1:end);    
    ifi=str2double(ifistr);
    
    tline=fgetl(f);
    nrepstr=tline(length('Repetitions: ')+1:end);    
    nrep=str2double(nrepstr); 
    
    tline=fgetl(f);
    tstr=tline(length('Speed: ')+1:end);    
    speed=str2double(tstr); 
        
    tline=fgetl(f);
    bw_str=tline(length('curtain_bw: ')+1:end);     
    curtain_bw=uint32(str2num(bw_str)); 

    tline=fgetl(f); %Baseline recording: 
    
    tline=fgetl(f);
    sf_str=tline(length('Spatial frequency: ')+1:end);     
    sp_freq=str2double(sf_str); 
    
    fclose(f);
end