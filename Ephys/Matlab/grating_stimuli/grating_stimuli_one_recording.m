% the code for analysing the recordings of grating stimulus 

% 07.01.2022
% O.Symonova
function one_rec_data = grating_stimuli_one_recording(logfile_fullname)
    [curfolder,~,~] = fileparts(logfile_fullname);        
    %find the closest pr file
    prfilename = find_pr_file_closest_date(logfile_fullname);
    one_rec_data={};
    if prfilename == -1
        warning(['Could not find the matching pr file for ',logfile_fullname]);
        return;
    end
    
    prfile_fullname=fullfile(curfolder,prfilename); 
    try
        one_rec_data = raw_ds_grating(prfile_fullname,logfile_fullname);
        close all;
    catch
        warning(['Something went wrong in analysis of ',logfile_fullname]);
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
   disp(min_et)
   if min_et>10
       disp('No matching pr file has been found');
       prfilename=-1;
   end      
end

function exp_data = raw_ds_grating(prfile_fullname,logfile_fullname)
% function  save_raw_ds_grating(rootDir, prfullname,logfullname, filedatestr) 
    %% get celltype and strain type from the path
    [folderpath,~,~] = fileparts(logfile_fullname);
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(folderpath);

    date_year=str2num(datestr(1:2));
    date_month=str2num(datestr(3:4));
     
    if date_year<21 || (date_year==21 && date_month<4)
        old_version=1; %in old version the frames_counter started from 1 hence the first red frame is on the 5th frame
    else
        old_version=0;
    end     

    %% open ephys data
    [Data, Text_header, filenameout, sampling_rate, minV, maxV]=openpr(prfile_fullname,0);
    
    %% list of params:   
    %channel in Data with red frames
    redchannel=2;  
    %value to threshold red frames
    red_threshold=4;
    %correct y axis scale to mV by a scale factor
    yscale=100;
    Data(:,1)=Data(:,1)*yscale;
        
    disp(['Analysing ',prfile_fullname]);
    [~,prname,~]=fileparts(prfile_fullname);   
    prfile=[prname,'.pr'];
    
    %add loading log file and saving parameters
    stimparam = load_grating_stimparam(logfile_fullname);
    
    %add processing multiple repetitions
    
    
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);%    
    startend=find_start_end_grating(Data, sampling_rate, old_version);
    bl_dur=stimparam.baseline_dur*sampling_rate;    
    %remove any red frame before the duration of the baseline, where no red
    %frame is supposed to be
    idk=find(startend(1,:)<bl_dur);
    startend(:,idk)=[];

    ndir=8;
    dir_vector = [0:45:359];
    baseline_duration=round(stimparam.baseline_dur/2*sampling_rate);
    nrep=size(startend,2)/ndir;
    
    min_grating_dur=min(startend(2,:)-startend(1,:));
    
    all_rec_n = size(startend,2);
    
    resps_all_raw  = zeros(nrep,ndir,min_grating_dur);
    baseline_raw  = zeros(nrep,ndir,baseline_duration);    
    for i=1:all_rec_n
        repi=floor((i-1)/ndir)+1;
        diri=i-(repi-1)*ndir;
        resps_all_raw(repi,diri,:)  = Data(startend(1,i):startend(1,i)+min_grating_dur-1,1);            
        baseline_raw(repi,diri,:)  = Data(startend(1,i)-baseline_duration+1:startend(1,i),1);        
    end
    
    
    av_baseline = mean(baseline_raw,3);
    resp_minus_baseliene=resps_all_raw;
    shifted_baseline=baseline_raw;  
    for ri=1:nrep
        resp_minus_baseliene(ri,:,:)=resp_minus_baseliene(ri,:,:)-av_baseline(ri,:);
        shifted_baseline(ri,:,:)=baseline_raw(ri,:,:)-av_baseline(ri,:);    
    end
    resp_minus_baseliene=reshape(resp_minus_baseliene,ndir,nrep,[]);
    shifted_baseline=reshape(shifted_baseline,ndir,nrep,[]);
    av_baseline = reshape(av_baseline,ndir,[]);
         
    
    exp_data.prname=prfile;
    exp_data.logname=logfile_fullname;
    exp_data.cellinfo.celltype=celltype;
    exp_data.cellinfo.cellgrp=cellgroup;    
    exp_data.cellinfo.strain_type =strain;
    exp_data.cellinfo.cellid =cellid;
    exp_data.cellinfo.datestr =datestr;        
    
    exp_data.siminfo.speed=stimparam.speed;
    exp_data.siminfo.sp_freq=stimparam.sp_freq;
    exp_data.siminfo.nrep=nrep;
    exp_data.siminfo.ndir=ndir;
    exp_data.siminfo.dir_vector=dir_vector;

    %raw data
    exp_data.all_repsp = resp_minus_baseliene;
    exp_data.all_baseline = shifted_baseline; 
    exp_data.mean_baseline = av_baseline; 
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
   drf4=median(drf(drf>0.016*sampling_rate));
   gap_min= 0.5*sampling_rate;
   
   %find big gaps in data (gap of 0.5sec is hardcoded in stimulus) 
   gap_rf_ind = find(drf>gap_min);
   lastind = [rft(gap_rf_ind); rft(end)];
   firstind= [rft(1); rft(gap_rf_ind+1)];
   if old_version_yn
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



function stimparam = load_grating_stimparam(logfullname)
% read the file with stimuli parameters
    f=fopen(logfullname);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'grating_stimuli'))
        disp('The log file does not contain information about grating stimuli.');
        fclose(f);
        return;
    end

    tline=fgetl(f); %start time of the recording    

    tline=fgetl(f); %end time of the recording
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparam.time=datevec(tstr,formatIn);
    
    tline=fgetl(f); %time per grating 
    tline=fgetl(f); %Flip rate: 
    
    tline=fgetl(f);
    tstr=tline(length('Repetitions: ')+1:end);    
    stimparam.nrep=str2double(tstr); 
    
    tline=fgetl(f);
    tstr=tline(length('Speed: ')+1:end);    
    stimparam.speed=str2double(tstr); 
        
    tline=fgetl(f);
    tstr=tline(length('Baseline recording: ')+1:end);     
    stimparam.baseline_dur=uint32(str2num(tstr)); 
    
    tline=fgetl(f);
    tstr=tline(length('Spatial frequency: ')+1:end);     
    stimparam.sp_freq=str2double(tstr);     
    
    fclose(f);
end
