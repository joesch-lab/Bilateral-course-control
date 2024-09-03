%the code for analysing the recordings during scanning rectangle stimulus
% 7.08.2020
% O.Symonova

parent_folder='C:\DATA\EPhys\fly_noise_test\CantonS_HSE\'
% parent_folder='C:\DATA\Data_for_Olga\fly_ephys\FlpND\201202';

% parent_folder='C:\DATA\Data_for_Olga\fly_ephys\FlpD\';
% parent_folder='C:\DATA\Data_for_Olga\fly_ephys\FlpD\200824\FlpD';
% parent_folder='C:\DATA\Data_for_Olga\fly_ephys';
% parent_folder='C:\DATA\Data_for_Olga\new_\';
% parent_folder='\\fs.ist.ac.at\drives\symonova\fs3-joeschgrp\Vika\EPhys\shakB_project\';
% parent_folder='\\fs.ist.ac.at\drives\symonova\fs3-joeschgrp\Vika\EPhys\shakB_project\';
% parent_folder='\\fs.ist.ac.at\drives\symonova\fs3-joeschgrp\Vika\EPhys\shakB_project\Repository - Copy\FlpD\VS\VS1_4\210302\scanning_rect_stimulus';
%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list
num_hor_bins=17;


%go thru the list, analize each log file with 'scanning_dot' in the title
for li =1:length(filelist)
    if  contains(filelist(li).name,'scanning','IgnoreCase', true) && contains(filelist(li).name,'.log','IgnoreCase', true)  
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        %find the closest pr file
        prfilename = find_pr_file_closest_date(logfile_fullname);
        folder=filelist(li).folder;
        prfile_fullname=fullfile(folder,filesep,prfilename);
        slashinds=strfind(logfile_fullname,filesep);
        dateind=find(slashinds>length(parent_folder),1,'first');
        lastind=min(length(folder),slashinds(dateind));
        filedatestr=folder(length(parent_folder)+1:lastind);
        filedatestr(filedatestr==filesep)=[];
        try
            compute_scanning_RF_OF_full_res(parent_folder, prfile_fullname,logfile_fullname, filedatestr, num_hor_bins);
        catch
            warning(['Something went wrong in analysis of ',logfile_fullname]);
        end
    end
end
      
 %function to find the pr file closes to the log file     
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


function compute_scanning_RF_OF_full_res(rootDir, prfullname,logfullname, filedatestr, num_hor_bins);
    %params for the analysis of the stimulus
    redchannel=2; 
    red_threshold =1; %value used to separate red frames from noisy background
    %in this interval there should only 1 red frame, others will be removed
    no_doubles_interval=0.02; %in seconds
    %offset to compute the responses of the cell to a frame
    resp_offset_sec = 0.0;
    % duration of one bin to compute the responses
    bin_duration_sec=0.08;
    

    %if filedatestr is empty take it as the creation date
    if isempty(filedatestr)
        file_info = dir(prfullname);
        filedatestr = datestr(file_info.date,'yyyymmdd');
    end
    [folderdata,prname,prext]=fileparts(prfullname);
    resname = [filedatestr,'_scanning_rect_',prname];     
    titlestr=strrep(resname,'_','\_');
    
    [filepath,name,ext] = fileparts(logfullname);
    resfolder=fullfile(filepath,'res');
    if ~exist(resfolder)
        mkdir(resfolder);
    end
    
    %% open ephys data
    %    [Data, Text_header, filenameout, sampling_rate]=openpr(prfullname,0);% openpr_flatten_translate(prfullname,0);
    [Data, Text_header, filenameout, sampling_rate]=openpr_flatten_translate(prfullname,1);
   
   %% find red frames, find the start and end of the stimulus, remove
    %%doubles and reconstruct missing frames
   
    Data = threshold_red_signal(Data, red_threshold);
    frame1=0.01*sampling_rate; %duration of one frame
    Data = remove_red_frames_repetions(Data, frame1);
    Data = insert_missing_red_frames(Data, sampling_rate,prname,resfolder);
  
   %get the red frames after the first cleanup
   rft=find(Data(:,2));  
   
   %find 1st and last frames of the stimulus, there are 3 consequtive red
   %frames around the stimulus presentation
   [rf1, rfe] = find_first_last_red_frames(Data, 3);
   rfe_time = rft(rfe);
     
   % remove repeatead red frames during stimulus
   no_doubles_interval_sr=no_doubles_interval*sampling_rate;
   Data = remove_red_frames_repetions(Data, no_doubles_interval_sr,rf1, rfe-1);
   %adjust the id of the last frame after removing duplicates
   rft=find(Data(:, redchannel));
   rfe=find(rft==rfe_time);
   
   ifired=median(rft(rf1+1:rfe)-rft(rf1:rfe-1));
   frame_duration=ifired/5;
   
      
   figure, plot(Data(:,2));
   hold on, plot([rft(rf1),rft(rf1)],[0,1],'g');
   hold on, plot([rft(rfe),rft(rfe)],[0,1],'r');
   
   
   %bkgResp=mean(Data(rfe1:rfee,1));%baseline at the end
   bkgResp=mean(Data(rft(rf1):rft(rfe),1)); %baseline as the mean of the recording
   Data(:,1)=Data(:,1)-bkgResp;
   
%% using lin interpolation find the timing of every frame
   nfr_act=(rfe-rf1+1)*5; %number of all frames of the stimulus
   allframesst=1:nfr_act;% %id of all frames
   rfst=1:5:nfr_act; %id of red frames
   rft_st=rft(rf1:rfe)'; %timing of red frames as recorded
   %liniar interpolation to find timing of all frames
   frametiming=uint32(round(interp1(rfst,rft_st,allframesst,'linear','extrap')));
   
    %% reconstruct frame array 
    [stim_arr, nrep]=reconstruct_scanning_rect_with_pause_fullres(logfullname);
    %stim_arr is [Nframes x 1 x 4] array, encodes position [x,y,dx,dy] of
    %the rectangle at each frame   
    
    %convolve full res picture with the voltage response then bin it as
    %needed, save full res in mat file
    
    %time intervals to compute voltage response to a frame
    nfrstim = size(stim_arr,1);    
    bintiming_offset=resp_offset_sec*sampling_rate;%0;
    bin_duration=bin_duration_sec*sampling_rate;
    num_frame_registered=min(nfrstim+1,length(frametiming));
    bin_edges_st=frametiming(1:num_frame_registered);
    bin_edges_en=frametiming(1:num_frame_registered)+bin_duration;
    bin_edges=[bin_edges_st;bin_edges_en]'+bintiming_offset;
   
    %voltage responses of a frame
    V_frame_resp= cell2mat(arrayfun(@(i) mean(Data(bin_edges(i,1):bin_edges(i,2),1)),1:num_frame_registered-1, 'UniformOutput', false)');
   
    %convolve each frame with the response
    RFx=zeros(342,608);
    RFy=zeros(342,608);
    RFx_pos=zeros(342,608);
    RFx_neg=zeros(342,608);
    RFy_pos=zeros(342,608);
    RFy_neg=zeros(342,608);
    
    for i=1:num_frame_registered-1
        [frdx,frdy] = make_frame(stim_arr(i,:),608,342);
        RFx=RFx+frdx.*V_frame_resp(i);
        RFy=RFy+frdy.*V_frame_resp(i);
        %compute individual components
        RFx_pos=RFx_pos+(frdx>0).*V_frame_resp(i);
        RFx_neg=RFx_neg+(frdx<0).*V_frame_resp(i);
        RFy_pos=RFy_pos+(frdy>0).*V_frame_resp(i);
        RFy_neg=RFy_neg+(frdy<0).*V_frame_resp(i);
    end
    RFx=RFx./nrep;
    RFy=RFy./nrep;
    
    RFx_pos=RFx_pos./nrep;
    RFx_neg=RFx_neg./nrep;
    RFy_pos=RFy_pos./nrep;
    RFy_neg=RFy_neg./nrep;
    
    [RFx_hres, RFy_hres]= upscale_scanning_dimensions(RFx,RFy);
    [RFx_pos_hres, RFy_pos_hres]= upscale_scanning_dimensions(RFx_pos,RFy_pos);
    [RFx_neg_hres, RFy_neg_hres]= upscale_scanning_dimensions(RFx_neg,RFy_neg);
    
    [bin_dx_av, bin_resp, bin_resp_deformed, bin_hor, bin_ver, bin_hor_sep, bin_ver_sep, el_bincenters, az_bincenters] = deform_plot_uniform_sampling(RFx_hres, RFy_hres, num_hor_bins, 'c', 0);
      
    matfile=fullfile(resfolder, [name,'.mat']);
    figbase = fullfile(resfolder,name); 
    save(matfile, 'V_frame_resp','bin_edges','Data', 'stim_arr', 'bin_dx_av','bin_resp','bin_resp_deformed','RFx_hres', 'RFy_hres', 'RFx_pos_hres', 'RFx_neg_hres', 'RFy_pos_hres', 'RFy_neg_hres', 'bin_hor','bin_ver','el_bincenters', 'az_bincenters', 'figbase','name', '-nocompression','-v7.3');
    makeRFfigure(bin_dx_av, bin_resp_deformed, figbase);
    makeRFfigure_with_histograms(bin_dx_av, bin_hor, bin_ver, figbase);
    makeRF_indiv_compfigure(RFx_pos_hres, RFx_neg_hres, RFy_pos_hres, RFy_neg_hres, num_hor_bins, figbase);
    close all;    
     
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
   
   if num_missed>1  %pehaps more than one red frame was missed
       [datafolder,resfolder]=fileparts(resfolder);
       prfullname=fullfile(datafolder,prfilename);
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




    



function [RFx_hres, RFy_hres]= upscale_scanning_dimensions(RFx,RFy)
    %find non-zero rows and columns
    sumrows=sum(RFx,2);
    sumcols=sum(RFx,1);
    nnzrows=find(sumrows~=0);
    nnzcols=find(sumcols~=0);
    RFxeffective=RFx(nnzrows,nnzcols);
    
    sumrows=sum(RFy,2);
    sumcols=sum(RFy,1);
    nnzrows=find(sumrows~=0);
    nnzcols=find(sumcols~=0);
    RFyeffective=RFy(nnzrows,nnzcols);    
    
    %upscale non-zero components
    maxrows=max(size(RFxeffective,1), size(RFyeffective,1));
    maxcols=max(size(RFxeffective,2), size(RFyeffective,2));
    
    RFx_hres=imresize(RFxeffective,[maxrows,maxcols],'nearest');
    RFy_hres=imresize(RFyeffective,[maxrows,maxcols],'nearest');    
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

function [st,en]=get_start_end_frames(frameseq)
    nfr=length(frameseq);
    st=zeros(nfr,1);
    en=zeros(nfr,1);
    seqnum=1;
    fi = frameseq(1);
    st(seqnum)=fi;
    for fii=2:nfr
        if frameseq(fii)~=fi+1 %new seq, need to finish the prev
            en(seqnum)=fi;
            seqnum=seqnum+1;        
            fi=frameseq(fii);
            st(seqnum)=fi;
        else
            fi=frameseq(fii);
        end    
    end
    en(seqnum)=frameseq(end);
    st(seqnum+1:end)=[];
    en(seqnum+1:end)=[];
end

function vargout=subplot_tight(m, n, p, margins, varargin)
    %% subplot_tight
    % A subplot function substitude with margins user tunabble parameter.
    %
    %% Syntax
    %  h=subplot_tight(m, n, p);
    %  h=subplot_tight(m, n, p, margins);
    %  h=subplot_tight(m, n, p, margins, subplotArgs...);
    %
    %% Description
    % Our goal is to grant the user the ability to define the margins between neighbouring
    %  subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the
    %  margins between subplots can reach 40% of figure area, which is pretty lavish. While at
    %  the begining the function was implememnted as wrapper function for Matlab function
    %  subplot, it was modified due to axes del;etion resulting from what Matlab subplot
    %  detected as overlapping. Therefore, the current implmenetation makes no use of Matlab
    %  subplot function, using axes instead. This can be problematic, as axis and subplot
    %  parameters are quie different. Set isWrapper to "True" to return to wrapper mode, which
    %  fully supports subplot format.
    %
    %% Input arguments (defaults exist):
    %   margins- two elements vector [vertical,horizontal] defining the margins between
    %        neighbouring axes. Default value is 0.04
    %
    %% Output arguments
    %   same as subplot- none, or axes handle according to function call.
    %
    %% Issues & Comments
    %  - Note that if additional elements are used in order to be passed to subplot, margins
    %     parameter must be defined. For default margins value use empty element- [].
    %  - 
    %
    %% Example
    % close all;
    % img=imread('peppers.png');
    % figSubplotH=figure('Name', 'subplot');
    % figSubplotTightH=figure('Name', 'subplot_tight');
    % nElems=17;
    % subplotRows=ceil(sqrt(nElems)-1);
    % subplotRows=max(1, subplotRows);
    % subplotCols=ceil(nElems/subplotRows);
    % for iElem=1:nElems
    %    figure(figSubplotH);
    %    subplot(subplotRows, subplotCols, iElem);
    %    imshow(img);
    %    figure(figSubplotTightH);
    %    subplot_tight(subplotRows, subplotCols, iElem, [0.0001]);
    %    imshow(img);
    % end
    %
    %% See also
    %  - subplot
    %
    %% Revision history
    % First version: Nikolay S. 2011-03-29.
    % Last update:   Nikolay S. 2012-05-24.
    %
    % *List of Changes:*
    % 2012-05-24
    %  Non wrapping mode (based on axes command) added, to deal with an issue of disappearing
    %     subplots occuring with massive axes.
    %% Default params
    isWrapper=false;
    if (nargin<4) || isempty(margins)
        margins=[0.04,0.04]; % default margins value- 4% of figure
    end
    if length(margins)==1
        margins(2)=margins;
    end
    %note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
    [subplot_col,subplot_row]=ind2sub([n,m],p);  
    height=(1-(m+1)*margins(1))/m; % single subplot height
    width=(1-(n+1)*margins(2))/n;  % single subplot width
    % note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
    subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
    subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   
    merged_height=subplot_rows*( height+margins(1) )- margins(1);   % merged subplot height
    merged_width= subplot_cols*( width +margins(2) )- margins(2);   % merged subplot width
    merged_bottom=(m-max(subplot_row))*(height+margins(1)) +margins(1); % merged subplot bottom position
    merged_left=min(subplot_col)*(width+margins(2))-width;              % merged subplot left position
    pos=[merged_left, merged_bottom, merged_width, merged_height];
    if isWrapper
       h=subplot(m, n, p, varargin{:}, 'Units', 'Normalized', 'Position', pos);
    else
       h=axes('Position', pos, varargin{:});
    end
    if nargout==1
       vargout=h;
    end
end