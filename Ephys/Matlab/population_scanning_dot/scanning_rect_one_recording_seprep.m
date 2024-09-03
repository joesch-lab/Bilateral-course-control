%the code for analysing the recordings during scanning rectangle stimulus
% 7.08.2020
% O.Symonova

function raw_data = scanning_rect_one_recording_seprep(logfile_fullname)
raw_data={};
%find the closest pr file
prfilename = find_pr_file_closest_date(logfile_fullname);
if prfilename == -1
    warning(['Could not find the matching pr file for ',logfile_fullname]);
    return;
end
[folder,~,~] = fileparts(logfile_fullname);
prfile_fullname=fullfile(folder,prfilename);
%     try
raw_data = compute_scanning_RF_OF_perfile_seprep(prfile_fullname,logfile_fullname);
%     catch
%         warning(['Something went wrong in analysis of ',logfile_fullname]);
%     end
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



function exp_data = compute_scanning_RF_OF_perfile_seprep(prfullname,logfullname)
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

%% remove outliers in voltage Data
%baseline had been already substracted when reading the file
Data =  remove_outliers(Data);

%% clean red channel
red_sygnal_raw=Data(:,2);
Data = clean_red_signal(Data, sampling_rate,prname,resfolder);

%params for the analysis of the stimulus
%offset to compute the responses of the cell to a frame
resp_offset_sec = 0.0;
% duration of one bin to compute the responses
bin_duration_sec=0.08;

%% using lin interpolation find the timing of every frame
n_redf=sum(Data(:,2)); %number of red frames is the sum of 1s in the cleaned red channel
nfr_act=n_redf*5+5; %number of all frames of the stimulus
allframesst=1:nfr_act;% %id of all frames
rft_st=find(Data(:,2)); %timing of red frames as recorded
rfst=1:5:length(rft_st)*5; %id of red frames
%linear interpolation to find timing of all frames
frametiming=uint32(round(interp1(rfst,rft_st,allframesst,'linear','extrap')));
ifi=median(frametiming(2:end)-frametiming(1:end-1))/sampling_rate;
%% reconstruct frame array
[stim_arr, nrep]=reconstruct_scanning_rect_with_pause_fullres(logfullname);
%stim_arr is [Nframes x 1 x 4] array, encodes position [x,y,dx,dy] of
%the rectangle at each frame

%time intervals to compute voltage response to a frame
%binning is done using the timing of stimulus frames
nfrstim = size(stim_arr,1);
nfrstim_all = nfrstim * nrep;
bintiming_offset=resp_offset_sec*sampling_rate;%0;
bin_duration=bin_duration_sec*sampling_rate;
num_frame_registered=min(nfrstim_all,length(frametiming));
bin_edges_st=frametiming(1:num_frame_registered);
bin_edges_en=frametiming(1:num_frame_registered)+bin_duration;
bin_edges=[bin_edges_st;bin_edges_en]'+bintiming_offset;

%voltage responses to a frame
%compute mean response during the frame presentation
%     V_frame_resp= cell2mat(arrayfun(@(i) mean(Data(bin_edges(i,1):bin_edges(i,2),1)),1:num_frame_registered, 'UniformOutput', false)')';
V_frame_resp= cell2mat(arrayfun(@(i) mean(Data(bin_edges(i,1):bin_edges(i,2),1)),1:num_frame_registered, 'UniformOutput', false)')';
%among all reps compute mean response to a frame
vrep=nan(nrep,nfrstim);
for ri=1:nrep
    frst=(ri-1)*nfrstim+1;
    fren=ri*nfrstim;
    if fren>num_frame_registered
        fr_en=num_frame_registered-(ri-1)*nfrstim;
        fren=num_frame_registered;
    else
        fr_en=nfrstim;
    end
    vrep(ri,1:fr_en)= V_frame_resp(frst:fren);
    if fr_en<nfrstim
        break;
    end
end
V_frame_mean = mean(vrep,'omitnan');

%within a rep find first and last frames of horiz and vert motion
f1_h=find(stim_arr(:,1,3)~=0,1,"first");
fn_h=find(stim_arr(:,1,3)~=0,1,"last");
f1_v=find(stim_arr(:,1,4)~=0,1,"first");
fn_v=find(stim_arr(:,1,4)~=0,1,"last");

%find the first frames of motion to the left, to the rigth
hmove=squeeze(stim_arr(f1_h:fn_h,1,3));
hleft=double(hmove>0);
f1_h_left=find(diff([0;hleft])>0);
f1_h_left=f1_h_left+f1_h-1;
hright=double(hmove<0);
f1_h_right=find(diff([0;hright])>0);
f1_h_right=f1_h_right+f1_h-1;

%find the first frames of motion up/down
vmove=squeeze(stim_arr(f1_v:fn_v,1,4));
vdown=double(vmove<0);
f1_v_down=find(diff([0;vdown])>0);
f1_v_down=f1_v_down+f1_v-1;
vup=double(vmove>0);
f1_v_up=find(diff([0;vup])>0);
f1_v_up=f1_v_up+f1_v-1;

params=read_scanning_rect_with_pause(logfullname);

%% print_red_frame_estimates
filereport=fullfile(resfolder,[prname,'_missing_frames_info.txt']);
t_elapsed = etime(params.tend,params.tstart); %elapsed seconds
nrf_etime=floor((t_elapsed-params.tbefore -params.tafter-6*ifi)*60/5);
nrf_frame_counter=floor(params.nfr/5);
if mod(params.nfr,5)==0
    nrf_frame_counter=nrf_frame_counter+1;
end
nrf_redsignal=length(rft_st);
f=fopen(filereport,'a');
fprintf(f,'RF from elapsed time: %d\n',nrf_etime);
fprintf(f,'RF from frame counter: %d\n',nrf_frame_counter);
fprintf(f,'RF from red signal: %d\n',nrf_redsignal);
fclose(f);




freeze_frames=double(params.freezefr);
%freezeing frames are before each motion to the left and down
freeze_start=[f1_h_left-freeze_frames;f1_v_down-freeze_frames];
%for each scanning line substract the baseline computed during the
%freezing frames
for ri=1:nrep
    frst=(ri-1)*nfrstim+1;
    for fri=1:length(freeze_start)
        freezi_start_t=frametiming(frst+freeze_start(fri)-1);
        freezi_end_t=frametiming(frst+freeze_start(fri)-1+freeze_frames);
        freezei_bg=mean(Data(freezi_start_t:freezi_end_t,1));

        %substract the baseline for the frames till the next freeze
        if fri==length(freeze_start)
            freezi_start_ii=nfrstim;
        else
            freezi_start_ii=freeze_start(fri+1)-1;
        end
        vrep(ri,freeze_start(fri):freezi_start_ii)=vrep(ri,freeze_start(fri):freezi_start_ii)-freezei_bg;
    end
end




%get the trace per rep
stop_flag=0;
for ri=1:nrep
    % save binned trace for horizontal and vertical sweep
    %convert frame index into time
    frst=(ri-1)*nfrstim+f1_h;
    t1=frametiming(frst);

    fren=(ri-1)*nfrstim+fn_h;
    if fren>num_frame_registered
        fren=num_frame_registered;
        stop_flag =1;
    end
    tn=frametiming(fren);
    resp_vals(ri).hor = Data(t1:tn,1);
    resp_vals(ri).hor_rf_clean = find(Data(t1:tn,2));
    resp_vals(ri).hor_rf = red_sygnal_raw(t1:tn);
    resp_vals(ri).hor_left_start = frametiming((ri-1)*nfrstim+f1_h_left)-t1+1;
    resp_vals(ri).hor_right_start = frametiming((ri-1)*nfrstim+f1_h_right)-t1+1;
    if stop_flag
        break;
    end

    %vertical trace
    frst=(ri-1)*nfrstim+f1_v;
    t1=frametiming(frst);

    fren=(ri-1)*nfrstim+fn_v;
    if fren>num_frame_registered
        fren=num_frame_registered;
        stop_flag =1;
    end
    tn=frametiming(fren);
    resp_vals(ri).ver = Data(t1:tn,1);
    resp_vals(ri).ver_rf_clean = find(Data(t1:tn,2));
    resp_vals(ri).ver_rf = red_sygnal_raw(t1:tn);
    resp_vals(ri).ver_down_start = frametiming((ri-1)*nfrstim+f1_v_down)-t1+1;
    resp_vals(ri).ver_up_start  = frametiming((ri-1)*nfrstim+f1_v_up)-t1+1;
    if stop_flag
        break;
    end
end

%get the scan lines of the stimulus
[scan_h_x, scan_h_y, scan_v_x, scan_v_y] = get_scanlines_from_stimarray(stim_arr);

nxh=length(scan_h_x);
nyh=length(scan_h_y);
nxv=length(scan_v_x);
nyv=length(scan_v_y);


%get full res RF
RFx=zeros(nyh,nxh,nrep);
RFy=zeros(nyv,nxv,nrep);
RFx_pos=zeros(nyh,nxh,nrep);
RFx_neg=zeros(nyh,nxh,nrep);
RFy_pos=zeros(nyv,nxv,nrep);
RFy_neg=zeros(nyv,nxv,nrep);
RFx_cellcount=zeros(nyh,nxh,nrep);
RFy_cellcount=zeros(nyv,nxv,nrep);
RFx_pos_cellcount=zeros(nyh,nxh,nrep);
RFx_neg_cellcount=zeros(nyh,nxh,nrep);
RFy_pos_cellcount=zeros(nyv,nxv,nrep);
RFy_neg_cellcount=zeros(nyv,nxv,nrep);

for i=1:nfrstim
    if stim_arr(i,1,3)~=0
        rowind=find(scan_h_y==stim_arr(i,1,2));
        colind=find(scan_h_x==stim_arr(i,1,1));
        for ri=1:nrep
            RFx(rowind,colind,ri)= RFx(rowind,colind,ri)+stim_arr(i,1,3).*vrep(ri,i);
            RFx_cellcount(rowind,colind,ri)=RFx_cellcount(rowind,colind,ri)+1;
            if stim_arr(i,1,3)>0
                RFx_pos(rowind,colind,ri)= RFx_pos(rowind,colind,ri) + vrep(ri,i);
                RFx_pos_cellcount(rowind,colind,ri)=RFx_pos_cellcount(rowind,colind,ri)+1;
            else
                RFx_neg(rowind,colind,ri)= RFx_neg(rowind,colind,ri) + vrep(ri,i);
                RFx_neg_cellcount(rowind,colind,ri)=RFx_neg_cellcount(rowind,colind,ri)+1;
            end
        end
    elseif stim_arr(i,1,4)~=0
        rowind=find(scan_v_y==stim_arr(i,1,2));
        colind=find(scan_v_x==stim_arr(i,1,1));
        for ri=1:nrep
            RFy(rowind,colind,ri)= RFy(rowind,colind,ri)+stim_arr(i,1,4).*vrep(ri,i);
            RFy_cellcount(rowind,colind,ri)=RFy_cellcount(rowind,colind,ri)+1;
            if stim_arr(i,1,4)>0
                RFy_pos(rowind,colind,ri)= RFy_pos(rowind,colind,ri) + vrep(ri,i);
                RFy_pos_cellcount(rowind,colind,ri)=RFy_pos_cellcount(rowind,colind,ri)+1;
            else
                RFy_neg(rowind,colind,ri)= RFy_neg(rowind,colind,ri) + vrep(ri,i);
                RFy_neg_cellcount(rowind,colind,ri)=RFy_neg_cellcount(rowind,colind,ri)+1;
            end
        end
    end
end

RFx=RFx./RFx_cellcount;
RFx_pos=RFx_pos./RFx_pos_cellcount;
RFx_neg=RFx_neg./RFx_neg_cellcount;
RFy=RFy./RFy_cellcount;
RFy_pos=RFy_pos./RFy_pos_cellcount;
RFy_neg=RFy_neg./RFy_neg_cellcount;

for ri=1:nrep
    resp_vals(ri).RFx=RFx(:,:,ri);
    resp_vals(ri).RFy=RFy(:,:,ri);
    resp_vals(ri).RFx_pos=RFx_pos(:,:,ri);
    resp_vals(ri).RFx_neg=RFx_neg(:,:,ri);
    resp_vals(ri).RFy_pos=RFy_pos(:,:,ri);
    resp_vals(ri).RFy_neg=RFy_neg(:,:,ri);

    %         [resp_vals(ri).RFx, resp_vals(ri).RFy] = upscale_scanning_dimensions_fullscreen(RFx(:,:,ri)/nxmove,RFy(:,:,ri)/nymove);
    %         [resp_vals(ri).RFx_pos, resp_vals(ri).RFy_pos] = upscale_scanning_dimensions_fullscreen(RFx_pos(:,:,ri)/nxpos,RFy_pos(:,:,ri)/nypos);
    %         [resp_vals(ri).RFx_neg, resp_vals(ri).RFy_neg] = upscale_scanning_dimensions_fullscreen(RFx_neg(:,:,ri)/nxneg,RFy_neg(:,:,ri)/nyneg);
    %
end


%   %% for each rep
%     % reconstruct RF
%     % save binned trace for horizontal and vertical sweep
%     % compute response distribution
%     nsamples=1000;
%     binframe=0.01667*sampling_rate;
%     resp_vals={};
%     minV=min(Data(:,1));
%     maxV=max(Data(:,1));
%     vedges=linspace(minV,maxV,100);
%
%
%     for ri=1:nrep
%         % save binned trace for horizontal and vertical sweep
%         resp_vals(ri).hor = vrep(ri,f1_h:fn_h);
%         resp_vals(ri).ver = vrep(ri,f1_v:fn_v);
%
%         %start and end of the rep
%         frst=(ri-1)*nfrstim+1;
%         fren=ri*nfrstim;
%         tst=round(frametiming(frst));
%         ten=round(frametiming(fren));
%         % compute response distribution
%         %trace of the whole rep
%         tracerep=Data(tst:ten,1);
%         npts=randi(ten-tst-binframe,1, nsamples); %rantom time pts to get responses
%         allresp = cell2mat(arrayfun(@(i) mean(Data(npts(i)+tst:npts(i)+tst+binframe,1)),1:nsamples, 'UniformOutput', false)')';
%         vcount = histcounts(allresp,vedges);
%         resp_vals(ri).resp_hist = vcount/nsamples;
%         resp_vals(ri).vedges = vedges;
%
%         %get full res RF
%         RFx=zeros(342,608);
%         RFy=zeros(342,608);
%         RFx_pos=zeros(342,608);
%         RFx_neg=zeros(342,608);
%         RFy_pos=zeros(342,608);
%         RFy_neg=zeros(342,608);
%
%         nxmove=0; nymove=0;
%         nxpos=0; nxneg=0;
%         nypos=0; nyneg=0;
%         for i=1:nfrstim
%             if stim_arr(i,1,3)~=0, nxmove=nxmove+1; end
%             if stim_arr(i,1,4)~=0, nymove=nymove+1; end
%             if stim_arr(i,1,3)>0, nxpos=nxpos+1; end
%             if stim_arr(i,1,3)<0, nxneg=nxneg+1; end
%             if stim_arr(i,1,4)>0, nypos=nypos+1; end
%             if stim_arr(i,1,4)<0, nyneg=nyneg+1; end
%
%             [frdx,frdy] = make_frame(stim_arr(i,:),608,342);
%             RFx(:,:)= RFx(:,:)+frdx.*vrep(ri,i);
%             RFy(:,:)= RFy(:,:)+frdy.*vrep(ri,i);
%             %compute individual components
%             RFx_pos(:,:)=RFx_pos(:,:)+(frdx>0).*vrep(ri,i);
%             RFx_neg(:,:)=RFx_neg(:,:)+(frdx<0).*vrep(ri,i);
%             RFy_pos(:,:)=RFy_pos(:,:)+(frdy>0).*vrep(ri,i);
%             RFy_neg(:,:)=RFy_neg(:,:)+(frdy<0).*vrep(ri,i);
%         end
%         [resp_vals(ri).RFx, resp_vals(ri).RFy] = upscale_scanning_dimensions_fullscreen(RFx/nxmove,RFy/nymove);
%         [resp_vals(ri).RFx_pos, resp_vals(ri).RFy_pos] = upscale_scanning_dimensions_fullscreen(RFx_pos/nxpos,RFy_pos/nypos);
%         [resp_vals(ri).RFx_neg, resp_vals(ri).RFy_neg] = upscale_scanning_dimensions_fullscreen(RFx_neg/nxneg,RFy_neg/nyneg);
%         resp_vals(ri).RFamp = vecnorm(cat(3,resp_vals(ri).RFx,resp_vals(ri).RFy),2,3);
%     end

%% save data
%file info
exp_data.prname=prname;
exp_data.logname=name;
exp_data.stim_info.nrep=nrep;
exp_data.stim_info.scan_h_x=scan_h_x;
exp_data.stim_info.scan_h_y=scan_h_y;
exp_data.stim_info.scan_v_x=scan_v_x;
exp_data.stim_info.scan_v_y=scan_v_y;

%raw data
exp_data.data.V_per_frame=V_frame_resp;
exp_data.data.resp_bins = bin_edges;
exp_data.data.v_reps=vrep;

%RF data
exp_data.RF = resp_vals;
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

function [scan_h_x, scan_h_y, scan_v_x, scan_v_y] = get_scanlines_from_stimarray(stim_arr)
f1_h=find(stim_arr(:,1,3)~=0,1,"first");
fn_h=find(stim_arr(:,1,3)~=0,1,"last");
f1_v=find(stim_arr(:,1,4)~=0,1,"first");
fn_v=find(stim_arr(:,1,4)~=0,1,"last");

hinds=find(stim_arr(f1_h:fn_h,1,3)~=0);
scan_h_x=unique(stim_arr(hinds,1,1))';
scan_h_y=unique(stim_arr(hinds,1,2))';

hinds=find(stim_arr(f1_v:fn_v,1,4)~=0);
hinds=hinds+f1_v-1;
scan_v_x=unique(stim_arr(hinds,1,1))';
scan_v_y=unique(stim_arr(hinds,1,2))';
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

 maxgap=max(drf);
 max_gap_seconds=maxgap/sampling_rate;

%fill the gap, insert the frame in the interval
n_inserted=0;
while ~isempty(missed_rf)
    for i=1:length(missed_rf)
        ti=rft(missed_rf(i));
        Data(ti+ifiredbefore,2)=1;
        n_inserted=n_inserted+1;
    end
    rft=find(Data(:,2)==1);
    drf=rft(2:end)-rft(1:end-1);
    missed_rf=find(drf>ifiredbefore*1.75);
end

%write a report about missing frames
if n_inserted>0
    filereport=fullfile(datafolder,resfolder,[prfilename,'_missing_frames_info.txt']);   
    f=fopen(filereport,'w');
    fprintf(f,'Number of gaps in red signal: %d\n',num_missed);
    fprintf(f,'Biggest gap btw red frames: %.4f\n',max_gap_seconds);
    fprintf(f,'Number of freame inserted: %d\n',n_inserted);
    fclose(f);
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
beforelast_frame=find(nbafter>=num_rf_around+1,1,'last');
%    [val,beforelast_frame]=max(nbafter(datalen2:end));
%    beforelast_frame=datalen2+ beforelast_frame -1;
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