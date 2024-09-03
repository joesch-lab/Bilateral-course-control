function exp_data=trace_on_of_flashes_one_withpreflashtrace_recording(logfile_fullname, flash_dur_s, tbefore_s, tafter_s, preflash_s)
        
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
    [Data, Text_header, filenameout, sampling_rate, minV, maxV]=openpr_flatten(prfile_fullname,0);
    
    %% list of params:   
    %correct y axis scale to mV by a scale factor
    yscale=100;
    Data(:,1)=Data(:,1)*yscale;
        
    disp(['Analysing ',prfile_fullname]);
    [~,prname,~]=fileparts(prfile_fullname);   
    prfile=[prname,'.pr'];
    
    %loading log file and saving parameters
    stimparam=load_ffflashes_stimparam(logfile_fullname);
    if stimparam.flash_dur < flash_dur_s || stimparam.tbefore<tbefore_s || stimparam.tafter<tafter_s
        return;
    end

    %threshold red channel to find real read frames
    redmax=max(Data(:,2));
    Data = threshold_red_signal(Data, redmax*0.75);%    
    startend=find_start_end_flash(Data, sampling_rate,prfile_fullname, old_version);
    nrep=size(startend,1);
    ondur=min(startend(:,2)-startend(:,1))+1;
    ondur=min(ondur,flash_dur_s*sampling_rate);
    if nrep==1 
        offdur=(tbefore_s+tafter_s)*sampling_rate;
    else
        offdur=min(startend(2:end,1)-startend(1:end-1,2))-1;
        offdur=min(offdur,(tbefore_s+tafter_s)*sampling_rate);
    end  

    flash_chunk=round(preflash_s*sampling_rate);
    flash_onset=flash_chunk+1;
    resps_all_on  = [];
    resps_all_off  = [];
    for i=1:nrep
        resps_all_on = cat(1,resps_all_on, Data(startend(i,1)-flash_chunk:startend(i,1)+flash_chunk-1,1)');        
        if i<nrep || nrep==1
            resps_all_off = cat(1,resps_all_off, Data(startend(i,2)+1-flash_chunk:startend(i,2)+1+flash_chunk-1,1)');
        end
    end

    exp_data.prname=prfile;
    exp_data.logname=logfile_fullname;
    exp_data.cellinfo.celltype=celltype;
    exp_data.cellinfo.cellgrp=cellgroup;    
    exp_data.cellinfo.strain_type =strain;
    exp_data.cellinfo.cellid =cellid;
    exp_data.cellinfo.datestr =datestr;        
        
    exp_data.siminfo.nrep=nrep;
    exp_data.siminfo.flash_color=stimparam.flash_color;
    exp_data.siminfo.bg_color=stimparam.bg_color;
    exp_data.siminfo.flash_chunk=flash_chunk;
    exp_data.siminfo.flash_chunk_s=preflash_s;
    exp_data.siminfo.flash_onset=flash_onset;   
        
    %raw data
    exp_data.trace.on_traces = resps_all_on;
    exp_data.trace.off_traces = resps_all_off;     
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

function startend=find_start_end_flash(Data, sampling_rate,prfile_fullname, old_version)
   %find the start and end time of flashes
   %red frames are shown only during off flasches
   
   %time of red frames
   rft=find(Data(:,2));     
   %find the time difference between two consecutive red frames
   drf=rft(2:end)-rft(1:end-1);
   
   %find big gaps in data 
   gap_min=0.1*sampling_rate;
   gap_rf_ind = find(drf>gap_min);

   %% if the stimulus was before 210421 red was leaking during the flashes
   % hence the flash is when there  is red,
   % in later stimulus red was shown during off-flashes
   if old_version
       firstind= [rft(1); rft(gap_rf_ind+1)];
       lastind = [rft(gap_rf_ind); rft(end)];
   else
       firstind= rft(gap_rf_ind)+1;
       lastind = rft(gap_rf_ind+1)-1;
   end

   startend=[firstind';lastind']';   

   figure, plot(Data(:,2));
   st=zeros(size(Data(:,2)));
   st(firstind)=0.5;
   hold on; plot(st,'LineWidth',2);
   en=zeros(size(Data(:,2)));
   en(lastind)=0.5;
   hold on; plot(en,'LineWidth',2);
   title(prfile_fullname,'Interpreter','none');
end

