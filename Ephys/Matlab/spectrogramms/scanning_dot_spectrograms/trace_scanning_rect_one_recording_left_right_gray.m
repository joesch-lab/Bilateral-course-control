%the code for extraction of the raw trace during the scanning dot stimulus
%when the dot was on the left side (not-stimulating) side of the screen

% 4.7.2022
% O.Symonova

% logfile_fullname='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\HS\HSN\220506\stimuli_scanning_rect_rescale_with_pause_2022_05_06_17_47_50.log';
% logfile_fullname = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpD\HS\HSN\220323\stimuli_scanning_rect_rescale_with_pause_2022_03_23_12_07_20.log';

function raw_trace = trace_scanning_rect_one_recording_left_right_gray(logfile_fullname) 
    raw_trace={};
    part_screen=1/3;
    %find the closest pr file
    prfilename = find_pr_file_closest_date(logfile_fullname);
    if prfilename == -1
        warning(['Could not find the matching pr file for ',logfile_fullname]);
        return;
    end
    [folder,~,~] = fileparts(logfile_fullname);           
    prfile_fullname=fullfile(folder,prfilename);                       
    try
        raw_trace = get_trace_left_right_part_screen(prfile_fullname,logfile_fullname,part_screen);        
    catch
        warning(['Something went wrong in analysis of ',logfile_fullname]);
    end        
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
               prfilename=currfile;
           end
       end
    end
    if min_et>10
       disp('No matching pr file has been found');
       prfilename=-1;
    end      
end



function exp_data = get_trace_left_right_part_screen(prfullname,logfullname, part_screen)
    disp(['Analysing ',prfullname]);
    %% open ephys data    
    [Data, Text_header, filenameout, sampling_rate]=openpr_flatten(prfullname,0);
    [~,prname,~]=fileparts(prfullname);
    
    %folder for the results
    [filepath,name,ext] = fileparts(logfullname);
    resfolder=fullfile(filepath,'res');
    if ~exist(resfolder)
        mkdir(resfolder);
    end
    
    %correct the amplitudes for specific recordings
    crazyVoltsSet=[[datetime(2022,02,15,11,05,0), datetime(2022,02,15,11,41,0)];...
                   [datetime(2022,01,17,00,00,0), datetime(2022,01,17,24,00,0)]];
    finfo=dir(prfullname);
    for checki=1:size(crazyVoltsSet,1)
        if finfo.date > crazyVoltsSet(checki,1) && finfo.date < crazyVoltsSet(checki,2) 
            Data(:,1)=Data(:,1)/5;
            break;
        end
    end
    
    %% clean red channel    
    Data = clean_red_signal(Data, sampling_rate,prname,resfolder);
    
    %% remove outliers in voltage Data
    Data =  remove_outliers(Data);
    bkgResp=median(Data(:,1)); %baseline as the median of the recording
    Data(:,1)=Data(:,1)-bkgResp;
   
    %% using lin interpolation find the timing of every frame    
    n_redf=sum(Data(:,2)); %number of red frames is the sum of 1s in the cleaned red channel
    nfr_act=n_redf*5+5; %number of all frames of the stimulus
    allframesst=1:nfr_act;% %id of all frames   
    rft_st=find(Data(:,2)); %timing of red frames as recorded
    rfst=1:5:length(rft_st)*5; %id of red frames
    %linear interpolation to find timing of all frames
    frametiming=uint32(round(interp1(rfst,rft_st,allframesst,'linear','extrap')));
   
    %% reconstruct frame array 
    [stim_arr, nrep]=reconstruct_scanning_rect_with_pause_fullres(logfullname);
    %stim_arr is [Nframes x 1 x 4] array, encodes position [x,y,dx,dy] of
    %the rectangle at each frame   

    %number of frames per rep from the reconstructed stimulus    
    nfrstim = size(stim_arr,1);       

    %find all frames when the rectangle was in the left side of the screen
    %find x value corresponding to the rectangle on the left
    allx=squeeze(stim_arr(:,1,1));
    minx=min(allx(:));
    maxx=max(allx(:));
    lastx=round(minx+(1-part_screen)*(maxx-minx));
    %find all frames when the rectangle was on the left
    idx=find(allx>lastx & (stim_arr(:,1,3)~=0 | stim_arr(:,1,4)~=0));
    %find the consequtive periods
    ddf=find(diff(idx)>1);
    fr_st_en_contra=[[idx(1),idx(ddf+1)'];[idx(ddf)',idx(end)]]';
    ntimes_h=size(fr_st_en_contra,1);
    fr_st_en_contra_h=[];
    fr_st_en_contra_v=[];
    for ti=1:ntimes_h
        id_st=fr_st_en_contra(ti,1);        
        if stim_arr(id_st,1,3)~=0 
            fr_st_en_contra_h=[fr_st_en_contra_h;fr_st_en_contra(ti,:)];
        end
        if stim_arr(id_st,1,4)~=0 
            fr_st_en_contra_v=[fr_st_en_contra_v;fr_st_en_contra(ti,:)];
        end
    end


    %find all frames when the rectangle was on the right
    lastx=round(minx+part_screen*(maxx-minx));    
    idx=find(allx<lastx & (stim_arr(:,1,3)~=0 | stim_arr(:,1,4)~=0));
    %find the consequtive periods
    ddf=find(diff(idx)>1);
    fr_st_en_ipsi=[[idx(1),idx(ddf+1)'];[idx(ddf)',idx(end)]]';
    ntimes_v=size(fr_st_en_ipsi,1);
    fr_st_en_ipsi_h=[];
    fr_st_en_ipsi_v=[];
    for ti=1:ntimes_v
        id_st=fr_st_en_ipsi(ti,1);
        if stim_arr(id_st,1,3)~=0    
            fr_st_en_ipsi_h=[fr_st_en_ipsi_h;fr_st_en_ipsi(ti,:)];
        end
        if stim_arr(id_st,1,4)~=0 
            fr_st_en_ipsi_v=[fr_st_en_ipsi_v;fr_st_en_ipsi(ti,:)];
        end
    end    
    
    %concatanate the trace when the rect on the left
    trace_all={};
    for ri=1:nrep
        trace_all(ri).reptrace_contra = cell2mat(arrayfun(@(i) Data(frametiming(fr_st_en_contra(i,1)):frametiming(fr_st_en_contra(i,2)),1)', 1:size(fr_st_en_contra,1),'UniformOutput',false));
        trace_all(ri).reptrace_contra_h = cell2mat(arrayfun(@(i) Data(frametiming(fr_st_en_contra_h(i,1)):frametiming(fr_st_en_contra_h(i,2)),1)', 1:size(fr_st_en_contra_h,1),'UniformOutput',false));
        trace_all(ri).reptrace_contra_v = cell2mat(arrayfun(@(i) Data(frametiming(fr_st_en_contra_v(i,1)):frametiming(fr_st_en_contra_v(i,2)),1)', 1:size(fr_st_en_contra_v,1),'UniformOutput',false));        
        trace_all(ri).reptrace_ipsi = cell2mat(arrayfun(@(i) Data(frametiming(fr_st_en_ipsi(i,1)):frametiming(fr_st_en_ipsi(i,2)),1)', 1:size(fr_st_en_ipsi,1),'UniformOutput',false));
        trace_all(ri).reptrace_ipsi_h = cell2mat(arrayfun(@(i) Data(frametiming(fr_st_en_ipsi_h(i,1)):frametiming(fr_st_en_ipsi_h(i,2)),1)', 1:size(fr_st_en_ipsi_h,1),'UniformOutput',false));
        trace_all(ri).reptrace_ipsi_v = cell2mat(arrayfun(@(i) Data(frametiming(fr_st_en_ipsi_v(i,1)):frametiming(fr_st_en_ipsi_v(i,2)),1)', 1:size(fr_st_en_ipsi_v,1),'UniformOutput',false));               
    end

    params=read_scanning_rect_with_pause(logfullname);
    tb_samples=params.tbefore*sampling_rate;
    ta_samples=params.tafter*sampling_rate;
    trace_grey= [Data(frametiming(1)-tb_samples:frametiming(1),1)', Data(frametiming(end):frametiming(end)+ta_samples,1)'];
    
    %% save data
    %file info
    exp_data.prname=prname;
    exp_data.logname=name;
    
    %raw data
    exp_data.data.trace=trace_all;
    exp_data.data.trace_gray = trace_grey;    
end


function Data = threshold_red_signal(Data, red_threshold)
    %Data is the Nx2 array, red frames are in the 2nd channel
    rft=find(Data(:,2)>red_threshold); %red frames
    Data(rft,2)=1; %set max values of red to one
    Data(Data(:,2)~=1,2)=0; % set to zero values smaller than 1
end

function Data = remove_red_frames_repetions(Data, rep_period, rf1, rfn)
    %if threre are red frames within time window rep_period after the 1st red
    %frame, they will be removed 
    rft=find(Data(:,2)); %red frames
     
    if ~exist('rf1','var')
        rf1=1;
    end
    if ~exist('rfn','var')
        rfn=length(rft);
    end
   
    for ri=rf1:rfn 
       Data(rft(ri)+1:rft(ri)+rep_period,2)=0;
    end   
end

function Data = insert_missing_red_frames(Data, sampling_rate, prfilename,resfolder)
 
   %get the red frames 
   rft=find(Data(:,2));     
   
   %find the time difference between two consecutive red frames
   drf=rft(2:end)-rft(1:end-1);
   
   %median intra-red-frame-interval
   ifiredbefore=median(drf);
   %find where there is a big gap between the frames
   missed_rf=find(drf>ifiredbefore*1.75);
   num_missed = length(missed_rf);
   [datafolder,resfolder]=fileparts(resfolder);
   prfullname=fullfile(datafolder,prfilename);
   
   if num_missed>1  %pehaps more than one red frame was missed       
       warning([prfullname,': ',num2str(num_missed), ' red frames were missed, will try to fix. Check the raw data.']);
   end
   
   %write a report about missing frames
   if num_missed>0
       filereport=fullfile(datafolder,resfolder,[prfilename,'_missing_frames_info.txt']);
       maxgap=max(drf);
       max_gap_seconds=maxgap/sampling_rate;
       f=fopen(filereport,'w');
       fprintf(f,'Number of the missed red frames: %d\n',num_missed);
       fprintf(f,'Biggest gap btw red frames: %.4f\n',max_gap_seconds);
       fclose(f);
   end
   
   %fill the gap, insert the frame in the interval
   while ~isempty(missed_rf)
        for i=1:length(missed_rf)
            ti=rft(missed_rf(i));
            Data(ti+ifiredbefore,2)=1;
        end
        rft=find(Data(:,2)==1);
        drf=rft(2:end)-rft(1:end-1);
        missed_rf=find(drf>ifiredbefore*1.75);
   end
end

function [first_red_frame_id, last_red_frame_id] = find_first_last_red_frames(Data, num_rf_around)      
   % find the beginning and end of the stimulus: num_rf_around consecutive red
   % frames mark the start and end, the stimulus starts with a redframe, 
   % after stimulus gray screen also starts with the red frame
   
   %get the red frames 
   rft=find(Data(:,2)); 
   ifi = median(rft(2:end)-rft(1:end-1))/5;
   %for each red frame find the number of red frames in a small 
   % neighborhood before and after that frame
   nbwin=int32((num_rf_around+1)*ifi);    
   datalen2=round(length(Data)/2);  
   %how many red frames in the neighborhood before/after each frame
   nbafter=movsum(Data(:,2),[0,nbwin]).*Data(:,2); 
   nbbefore = movsum(Data(:,2),[nbwin,0]).*Data(:,2);
   
   %find the max in the fist half of the recording
   [val,first_frame_time]=max(nbbefore(1:datalen2));
   first_red_frame_id=find(rft==first_frame_time);
   %find the max in the second half of the recording
   [val,beforelast_frame]=max(nbafter(datalen2:end));
   beforelast_frame=datalen2+ beforelast_frame -1;
   %find the last redframe smaller than beforelast_frame
   last_red_frame_id=find(rft<beforelast_frame,1,'last');
      
%    figure, plot(Data(:,2));
%    hold on; plot([rft(first_red_frame_id),rft(first_red_frame_id)],[0,1],'g');
%    hold on; plot([rft(last_red_frame_id),rft(last_red_frame_id)],[0,1],'r');   
   
end




   

function [frdx,frdy] = make_frame(stim_arr_i,screen_xres, screen_yres)
    %convert a condensed frame representation [ndots x [xpos,ypos,da, dx]]
    % to the full frame representation 
    frdx=zeros(screen_yres,screen_xres);
    frdy=zeros(screen_yres,screen_xres);
    npts=size(stim_arr_i,1);
    for i=1:npts
        frdx(stim_arr_i(i,2),stim_arr_i(i,1))=stim_arr_i(i,3);
        frdy(stim_arr_i(i,2),stim_arr_i(i,1))=stim_arr_i(i,4);
    end    
end

function clean_data = clean_red_signal(Data, sampling_rate,prname,resfolder)
   %params for the analysis of the stimulus    
    red_threshold =1; %value used to separate red frames from noisy background
    %in this interval there should only 1 red frame, others will be removed
    no_doubles_interval=0.02; %in seconds
   
   %% find red frames, find the start and end of the stimulus, remove
    %%doubles and reconstruct missing frames
   
    clean_data = threshold_red_signal(Data, red_threshold);
    frame1=0.01*sampling_rate; %duration of one frame
    clean_data = remove_red_frames_repetions(clean_data, frame1);
    clean_data = insert_missing_red_frames(clean_data, sampling_rate,prname,resfolder);
  
   %get the red frames after the first cleanup
   rft=find(clean_data(:,2));  
   
   %find 1st and last frames of the stimulus, there are 3 consequtive red
   %frames around the stimulus presentation
   [rf1, rfe] = find_first_last_red_frames(clean_data, 3);  
   rfe_time = rft(rfe);
   rf1_time = rft(rf1);
     
   % remove repeatead red frames during stimulus
   no_doubles_interval_sr=no_doubles_interval*sampling_rate;
   clean_data = remove_red_frames_repetions(clean_data, no_doubles_interval_sr,rf1, rfe-1); 
   %remove any red signal before the first red frame
   clean_data(1:rf1_time-1,2) = 0;
   %remove any red signal after the last red frame
   clean_data(rfe_time+1:end,2) = 0;
end