% the code for analysing the recordings of grating stimulus 

% 07.01.2022
% O.Symonova
function save_grating_raw_data(logfile_fullname)  
    %find the closest pr file
    prfilename = find_pr_file_closest_date(logfile_fullname);
    if prfilename == -1
        warning(['Could not find the matching pr file for ',logfile_fullname]);
        return;
    end
    [folder,logname,ext] = fileparts(logfile_fullname);           
    prfile_fullname=fullfile(folder,prfilename);            
    [~,filedatestr,~]  = fileparts(folder);            
    try
        save_raw_ds_grating(folder, prfile_fullname,logfile_fullname, filedatestr);
        close all;
    catch
        warning(['Something went wrong in analysis of ',logfile_fullname]);
    end        
end
   
function  save_raw_ds_grating(rootDir, prfullname,logfullname, filedatestr) 
    %% get celltype and strain type from the path
     [filepath,name,ext] = fileparts(logfullname);
      [filepath_1,datestr,~] = fileparts(filepath);
     [filepath_2,cell_type,~] = fileparts(filepath_1);
     [filepath_3,cell_grp,~] = fileparts(filepath_2);
     [filepath_4,strain_type,~] = fileparts(filepath_3);
     
     date_year=str2double(datestr(1:2));
     date_month=str2double(datestr(3:4));
     
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
    %correct y axis scale to mV by a scale factor
    yscale=100;
    Data(:,1)=Data(:,1)*yscale;
    
    %% init output folder and file names
    %if filedatestr is empty take it as the creation date
    if isempty(filedatestr)
        file_info = dir(prfullname);
        filedatestr = datestr(file_info.date,'yyyymmdd');
    end
    
    disp(['Analysing ',prfullname]);
    [~,prname,~]=fileparts(prfullname);
       
    resfolder=fullfile(filepath,'res');
    if ~exist(resfolder)
        mkdir(resfolder);
    end
    
    %add loading log file and saving parameters
    stimparam = load_grating_stimparam(logfullname);
    %add processing multiple repetitions
    
    
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);%    
    startend=find_start_end_grating(Data, sampling_rate, old_version);
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
    resp_minus_baseliene=resps_all_raw-av_baseline;
    shifted_baseline=baseline_raw-av_baseline;
    av_raw_responses=mean(resp_minus_baseliene,1);
    av_baseline=mean(shifted_baseline,1);
    std_raw_responses=std(resp_minus_baseliene,1);
    std_baseline=std(shifted_baseline,1);
    
    exp_data.prname=prname;
    exp_data.logname=name;
    exp_data.cell_type=cell_type;
    exp_data.cell_grp=cell_grp;    
    exp_data.strain_type =strain_type;
    
    exp_data.timestamp=stimparam.time;
    exp_data.speed=stimparam.speed;
    exp_data.sp_freq=stimparam.sp_freq;
    exp_data.nrep=nrep;
    exp_data.ndir=ndir;
    exp_data.dir_vector=dir_vector;
    %raw data
    exp_data.all_repsp = resps_all_raw;
    exp_data.all_baseline= baseline_raw;    
    %average across reps
    exp_data.av_responses=av_raw_responses;
    exp_data.av_baseline=av_baseline;
    exp_data.std_raw_responses=std_raw_responses;
    exp_data.std_baseline=std_baseline;
    exp_data.shifted_baseline=shifted_baseline;
    
    matfile=fullfile(resfolder, [name,'_raw_ds.mat']);
    %save the data structure
    save(matfile, 'exp_data', '-nocompression','-v7.3');    
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
