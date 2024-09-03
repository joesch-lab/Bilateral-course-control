% the code for analysing the recordings of velocity tuning grating stimulus

% 07.02.2022
% O.Symonova
function raw_data = vel_tuning_one_recording(logfile_fullname) 
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
        raw_data = get_raw_vel_tuning(prfile_fullname,logfile_fullname);        
    catch
        warning(['Something went wrong in analysis of ',logfile_fullname]);
    end        
end
   

function  exp_data = get_raw_vel_tuning(prfullname,logfullname) 
    %% open ephys data
    [Data, ~, ~, sampling_rate, ~, ~]=openpr(prfullname,0);
    
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
  
    %value to threshold red frames
    red_threshold=1;
  
    disp(['Analysing ',prfullname]);
    [~,prname,~]=fileparts(prfullname);
    
    [filepath,name,~] = fileparts(logfullname);
    resfolder=fullfile(filepath,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
        
    %threshold red channel to find real read frames
    Data = threshold_red_signal(Data, red_threshold);
    
    %get param from the log file
    stimparam = read_velocity_tuning_stimuli_log(logfullname);
    
    %find start and end frames of moving curtain
    startend=find_start_end_moving_curtain(Data, sampling_rate);
    %duration of the baseline at least 1 sec before the motion
    bl_dur=1*sampling_rate;    
    %remove any red frame before the duration of the baseline, where no red
    %frame is supposed to be
    idk=find(startend(1,:)<bl_dur);
    startend(:,idk)=[];
     
    
    %remove spurious red frames 
    tmin=0.25*stimparam.tmove*sampling_rate;
    indsremove=[];
    for ii=1:size(startend,2)
        if startend(2,ii)-startend(1,ii)<tmin
            indsremove=[indsremove,ii];
        end
    end
    startend(:,indsremove)=[];

    nrep=stimparam.nrep;
    
    %traces between start and end frames
    resp_dur=min(startend(2,:)-startend(1,:));
    speed_estimate=floor(size(startend,2)/2/nrep);
    nspeeds=min(speed_estimate, length(stimparam.speed_array));
    nruns=nspeeds*2*nrep;
    
    %trim potential red contamination
    startend=startend(:,1:nruns);
    
    
    %for each grating 2sec static, t sec moving, 1 sec static
    % take 1.5 sec of the baseline before each grating
    bl_duration=round(1*sampling_rate);
    ncondintions=nspeeds*2;
    all_traces_info={};
    mean_traces=zeros(nrep,ncondintions,resp_dur);   
    minv=inf;
    maxv=-inf;
    for i=1:nruns
        %mean baseline right before the motion
        start_bl_ind=max(1,startend(1,i)-bl_duration);
        baseline_i=mean(Data(start_bl_ind:startend(1,i)-1,1));
        all_traces_info(i).baseline = baseline_i;
        all_traces_info(i).trace=[Data(startend(1,i):startend(1,i)+resp_dur-1,1)-baseline_i]';
        all_traces_info(i).baseline_trace=[Data(start_bl_ind:startend(1,i)-1,1)]';
        minv=min(minv,min(all_traces_info(i).trace));
        maxv=max(maxv,max(all_traces_info(i).trace));
        ii=ceil(i/ncondintions);
        jj=i-(ii-1)*ncondintions;
        %(reps x conditions x time)       
        if stimparam.isHS 
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
        all_traces_info(i).speed = stimparam.speed_array(stimparam.speed_indices(ceil(i/2)));        
    end
    
%     %mean responses
%     mean_traces_info={};    
%     for i=1:ncondintions
%         %get all reps with the same speed and dir values        
%         mean_traces_info(i).mean_trace=squeeze(mean(mean_traces(:,i,:),1)); %mean trace across reps
%         mean_traces_info(i).std_trace=squeeze(std(mean_traces(:,i,:),1)); %std across reps
%         mean_traces_info(i).mean_val=mean(mean_traces_info(i).mean_trace);
%         mean_traces_info(i).speed=all_traces_info(i).speed;
%         mean_traces_info(i).dir=all_traces_info(i).dir;
%         mean_traces_info(i).nsamples=nrep;
%     end
              
    exp_data.prname=prname;
    exp_data.logname=name;
%     exp_data.folder=filepath;
%     
%     [strain,cellgroup,celltype,date_str] = cell_info_from_path(filepath);    
%     exp_data.cell_type=celltype;
%     exp_data.cell_grp=cellgroup;    
%     exp_data.strain_type =strain;
%     exp_data.timestamp=date_str;    
    
    exp_data.tmove=stimparam.tmove;
    exp_data.sp_freq=stimparam.sp_freq;
    exp_data.nrep=nrep;
    exp_data.isHS=stimparam.isHS;

    %data
    exp_data.raw_data = all_traces_info;
%     exp_data.mean_vals = mean_traces_info; 
    
%     matfile=fullfile(resfolder, [name,'_res_vel_tuning.mat']);
%     %save the data structure
%     save(matfile, 'exp_data', '-nocompression','-v7.3');  
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

