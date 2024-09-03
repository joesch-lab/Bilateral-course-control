function exp_data=trace_pd_nd_grating_one_recording(logfile_fullname)
% function  save_raw_ds_grating(rootDir, prfullname,logfullname, filedatestr) 
    %% get celltype and strain type from the path
    exp_data={};
    [folderpath,~,~] = fileparts(logfile_fullname);
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(folderpath);

    date_year=str2double(datestr(1:2));
    date_month=str2double(datestr(3:4));
     
    if date_year<21 || (date_year==21 && date_month<4)
        old_version=1; %in old version the frames_counter started from 1 hence the first red frame is on the 5th frame
    else
        old_version=0;
    end     

    %% get the matching pr file
    prfile_fullname = find_pr_file_closest_date(logfile_fullname);
    if prfile_fullname==-1
        return;
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
     
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);%    
    startend=find_start_end_grating(Data, sampling_rate, old_version);
    bl_dur=stimparam.baseline_dur*sampling_rate;    
    %remove any red frame before the duration of the baseline, where no red
    %frame is supposed to be
    idk=find(startend(1,:)<bl_dur);
    startend(:,idk)=[];

    if isempty(stimparam.VSHS)
        ndir=8;
        dir_vector = [0:45:359];
        if contains(cellgroup,'HS')
            pd_idx=[1];
            nd_idx=[5];
        else %VS cell
            pd_idx=[3];
            nd_idx=[7];
        end
    else
        ndir=2;
        if contains(stimparam.VSHS,'HS')
            if ~contains(cellgroup,'HS')
                return; %skip HS files for non HS cells
            end
            dir_vector = [0,180];
            pd_idx=[1];
            nd_idx=[2];
        else %VS cell
            if ~contains(cellgroup,'VS')
                return; %skip VS files for non VS cells
            end
            dir_vector = [90,270];
            pd_idx=[1];
            nd_idx=[2];
        end
    end

    baseline_duration=min(1,single(stimparam.baseline_dur))*sampling_rate;
    nrep=size(startend,2)/ndir;
    
    min_grating_dur=min(startend(2,:)-startend(1,:));
    
    all_rec_n = size(startend,2);
    
%     resps_all_pd  = zeros(nrep*length(pd_idx),min_grating_dur);
%     resps_all_nd  = zeros(nrep*length(nd_idx),min_grating_dur);
%     baseline_pd  = zeros(nrep*length(pd_idx),baseline_duration);    
%     baseline_nd  = zeros(nrep*length(nd_idx),baseline_duration);   
    resps_all_pd  = [];
    resps_all_nd  = [];
    baseline_pd  = [];
    baseline_nd  = [];
    for i=1:all_rec_n
        repi=floor((i-1)/ndir)+1;
        diri=i-(repi-1)*ndir;
        if ismember(diri,pd_idx)
            resps_all_pd = cat(1,resps_all_pd, Data(startend(1,i):startend(1,i)+min_grating_dur-1,1)');
            baseline_pd = cat(1,baseline_pd, Data(startend(2,i):startend(2,i)+baseline_duration,1)');
        elseif ismember(diri,nd_idx)
            resps_all_nd = cat(1,resps_all_nd, Data(startend(1,i):startend(1,i)+min_grating_dur-1,1)');             
            baseline_nd = cat(1,baseline_nd, Data(startend(2,i):startend(2,i)+baseline_duration,1)');
        end
    end          
    
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
    exp_data.siminfo.pd_ind=pd_idx;
    exp_data.siminfo.nd_ind=nd_idx;

    %raw data
    exp_data.trace.pd_traces = resps_all_pd;
    exp_data.trace.pd_baseline = baseline_pd; 
    exp_data.trace.nd_traces = resps_all_nd;
    exp_data.trace.nd_baseline = baseline_nd;     
end



 %function to find the pr file closes to the log file     
function prfilename = find_pr_file_closest_date(logfilename)
    [folder,~,~]=fileparts(logfilename);    
    %date info from log file
    log_info = dir(logfilename);
    log_datevec=datevec(log_info.date);
    D = dir(folder);
    min_et=3600;
    prfilename=-1;
    %find the closest pr file by time
    for k = 3:length(D) % avoid using the first ones
       currfile = D(k).name; %file name
       if contains(currfile,'.pr','IgnoreCase', true)
           et=abs(etime(datevec(D(k).date),log_datevec));
           if et<min_et
               min_et=et;
               prfilename=fullfile(D(k).folder,currfile);
           end
       end
    end
    if min_et>10
       disp('No matching pr file has been found');
       prfilename=-1;
    end      
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

