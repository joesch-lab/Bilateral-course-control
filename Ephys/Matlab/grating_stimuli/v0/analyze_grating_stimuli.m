% the code for analysing the recordings of grating stimulus of different
% contrasts

% 02.03.2021
% O.Symonova

parent_folder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND DB331\HS\HSE\220209\';
% parent_folder='\\istsmb3.ist.local\joeschgrp\Vika\EPhys\shakB_project\Repository - Copy\FlpD\HS\HSS\200205\curtain_stimulus\';
% C:\DATA\Data_for_Olga\fly_ephys\FlpND\201202\copy\curtain_stimulus

%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list


%go thru the list, analize each log file with 'scanning_dot' in the title
for li =1:length(filelist)
    if  contains(filelist(li).name,'grating_stimuli','IgnoreCase', true) && contains(filelist(li).name,'.log','IgnoreCase', true)  
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
            compute_grating(parent_folder, prfile_fullname,logfile_fullname, filedatestr);
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


function  compute_grating(rootDir, prfullname,logfullname, filedatestr) 
     %% get celltype and strain type from the path
     [filepath,name,ext] = fileparts(logfullname);
     [filepath_1,thisfiledate,~] = fileparts(filepath);
     [filepath_2,cell_type,~] = fileparts(filepath_1);
     [filepath_3,cell_grp,~] = fileparts(filepath_2);
     [filepath_4,strain_type,~] = fileparts(filepath_3);
     
     date_year=str2double(thisfiledate(1:2));
     date_month=str2double(thisfiledate(3:4));
     
     if date_year<21 || (date_year==21 && date_month<4)
         old_version=1; %in old version the frames_counter started from 1 hence the first red frame is on the 5th frame
     else
         old_version=0;
     end
     
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
    
    %add loading log file and saving parameters
    %add processing multiple repetitions
    
    
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);%    
    startend=find_start_end_grating(Data, sampling_rate, old_version);
    all_rec_n = size(startend,2);

    resps_all  = zeros(all_rec_n,1);
    for i=1:all_rec_n
        resps_all(i)  = mean(Data(startend(1,i):startend(2,i),1));
    end
    baseline_duration=1*sampling_rate;
    baseline_resps  = zeros(all_rec_n,1);
    for i=1:all_rec_n
        baseline_resps(i)  = mean(Data(startend(1,i)-baseline_duration:startend(1,i)-1,1));
    end
    resps_minus_baseline = resps_all-baseline_resps;
    figbase =fullfile(resfolder, [name,'_fig_']);
    make_polar_plot_direction_selectivity(av_total, mean_traces, baseline_activity, direction_vector, title_str, figbase)
  
    
    matfile=fullfile(resfolder, [name,'_res.mat']);

    %save variables
    save(matfile, 'resps_all', 'baseline_resps','resps_minus_baseline', '-nocompression','-v7.3');    
end


function Data = threshold_red_signal(Data, red_threshold)
    %Data is the Nx2 array, red frames are in the 2nd channel
    rft=find(Data(:,2)>red_threshold); %red frames
    Data(rft,2)=1; %set max values of red to one
    Data(Data(:,2)~=1,2)=0; % set to zero values smaller than 1
end

function startend=find_start_end_grating(Data, sampling_rate, old_version_yn)
   %find the start and end time of grating
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
   if old_version_yn
       %find the gap between every 5th frame
       rf4ins = find(drf>0.01667*2*sampling_rate);
       drf4=median(drf(rf4ins));   
       firstind= firstind-drf4;
   end
   
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
    thetaticks(direction_vector);
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