function raw_data = contrast_levels_one_recording(logfile_fullname) 
% the code for analysing the recordings of grating stimulus of different
% contrasts, returns raw traces during grating presentation and baseline
% values

% 18.02.2022
% O.Symonova

    raw_data={};
    %find the closest pr file
    prfilename = find_pr_file_closest_date(logfile_fullname);
    if prfilename == -1
        warning(['Could not find the matching pr file for ',logfile_fullname]);
        return;
    end
    [folder,~,~] = fileparts(logfile_fullname);           
    prfile_fullname=fullfile(folder,prfilename);                       
    try
        raw_data = compute_grating_contrasts_response(prfile_fullname,logfile_fullname);        
    catch
        warning(['Something went wrong in analysis of ',logfile_fullname]);
    end        
end


function raw_data =  compute_grating_contrasts_response(prfullname,logfullname) 
    raw_data = {};
    %% open ephys data
    [Data, ~, ~, sampling_rate, ~ , ~]=openpr(prfullname,1);  
    slashinds=strfind(prfullname,filesep);
    title(prfullname(slashinds(end-3):end));

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

    %% list of params:   
    %channel in Data with red frames
    redchannel=2;  
    %value to threshold red frames
    red_threshold=4;
        
    disp(['Analysing ',prfullname]);
    [~,prname,~]=fileparts(prfullname);
        
    
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);
%     minresp=min(Data(:,1));
%     Data(:,1)=Data(:,1)-minresp;
    
    stimparam = load(logfullname);
    nconfigs = size(stimparam.allconfigs,1);
    nrep = stimparam.repetitions;
           
    %find start and end frames of moving curtain
    baseline_dur = 1*sampling_rate;
    startend=find_start_end_grating(Data, sampling_rate);
    %remove any red frame before the duration of the baseline, where no red
    %frame is supposed to be
    idk=find(startend(1,:)<baseline_dur);
    startend(:,idk)=[];
    all_rec_n = size(startend,2);

%     %%check the extraction of the start and end of the grating    
%     figure, plot(Data(:,1)); hold on; plot(Data(:,2)); hold on;
%     for i=1:size(startend,2)
%         plot([startend(1,i), startend(1,i)], [0,1], 'g'); hold on;
%         plot([startend(2,i), startend(2,i)], [0,1], 'k'); hold on;
%     end


    mindur = min(startend(2,:)-startend(1,:));
    traces = zeros(all_rec_n,mindur);
    baseline_traces = zeros(all_rec_n, 5000);
    baseline = zeros(all_rec_n,1);
    for i=1:all_rec_n
        start_bl=max(1,startend(1,i)-baseline_dur);
        baseline(i)=mean(Data(start_bl:startend(1,i)-1,1));
        traces(i,:)=Data(startend(1,i):startend(1,i)+mindur-1,1)-baseline(i)';
        baseline_traces(i,:)=Data(startend(1,i)-5000:startend(1,i)-1,1)';
    end

    %every 1:nconfigs:end runs of the stimulus is the gray screen remove it
    %from traces
    if stimparam.withgray_yn
        warning('Analyzing conrast grating stimulus with gray screen. Check the values');
        gr_ind=1:nconfigs:all_rec_n;
        non_gray = setdiff(1:all_rec_n,gr_ind);
        baseline = baseline(non_gray);
        traces=traces(non_gray,:);
    end

    raw_data.baseline=baseline;
    raw_data.trace=traces;
    raw_data.baseline_traces=baseline_traces;
    [~,prname,prext]=fileparts(prfullname);
    [~,logname,logext]=fileparts(logfullname);
    raw_data.prname=[prname,prext]; 
    raw_data.logname=[logname,logext]; 
end

function Data = threshold_red_signal(Data, red_threshold)
    %Data is the Nx2 array, red frames are in the 2nd channel
    rft=find(Data(:,2)>red_threshold); %red frames
    Data(rft,2)=1; %set max values of red to one
    Data(Data(:,2)~=1,2)=0; % set to zero values smaller than 1
end

function startend=find_start_end_grating(Data, sampling_rate)
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
   
   startend=[firstind';lastind'];
   
%    figure, plot(Data);
%    plot_pair_frames(firstind,lastind);   
end


% function plot_pair_frames(timest,timeen)
%     %assume the figure with the ephys activity is the current
%     hold on 
%     plot([timest, timest],[0,1], 'm'); %start frame
%     hold on
%     plot([timeen, timeen],[0,1], 'k'); %end
%     hold off
% end

