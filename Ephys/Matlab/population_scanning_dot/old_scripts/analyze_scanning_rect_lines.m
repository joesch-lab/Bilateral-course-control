%the code for analysing the recordings during scanning rectangle stimulus
% 7.08.2020
% O.Symonova

parent_folder='C:\Users\User\Desktop\scanning_rect_random_dots\';
parent_folder='D:\Repository\FlpD\HS\HSE\200824';

% parent_folder='C:\DATA\Data_for_Olga\new_\';

%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

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
            compute_scanning_RF_OF(parent_folder, prfile_fullname,logfile_fullname, filedatestr);
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


function compute_scanning_RF_OF(rootDir, prfullname,logfullname, filedatestr)   
    %if filedatestr is empty take it as the creation date
    if isempty(filedatestr)
        file_info = dir(prfullname);
        filedatestr = datestr(file_info.date,'yyyymmdd');
    end
    slashinds=strfind(prfullname,filesep);
    dotind=strfind(prfullname,'.');
    prname=prfullname(slashinds(end)+1:dotind(end)-1);
    %name of the files with the results
    resname = [filedatestr,'_scanning_rect_',prname];     
    titlestr=strrep(resname,'_','\_');
%% open ephys data
%    [Data, Text_header, filenameout, fs]=openpr(prfullname,0);% openpr_flatten_translate(prfullname,0);
    [Data, Text_header, filenameout, fs]=openpr_flatten_translate(prfullname,0);

%    [Data, Text_header, filenameout, fs, minV, maxV]=openpr_recent_noise(prfullname,1);
   Data(1:15000,:)=[];
   redchannel=2;     

   %% find red frames, find the start and end of the stimulus, remove
%%doubles and reconstruct missing frames
   
   rft=find(Data(:,2)>1); %red frames
   frame1=0.01*fs; %duration of one frame
   Data(rft,2)=1; %set max values of red to one
   Data(Data(:,2)~=1,2)=0; % set to zero values smaller than 1
   
   %remove redframes that have been registered multiple times,
   % i.e. within the duration of one frame;
   for ri=1:length(rft) 
       Data(rft(ri)+1:rft(ri)+frame1,2)=0;
   end
   
   %get the red frames after the first cleanup
   rft=find(Data(:,2));  
   
   %reconstruct missing red frames
   %find the time difference between two consecutive red frames
   drf=rft(2:end)-rft(1:end-1);
   %median intra-red-frame-interval
   ifiredbefore=median(drf);
   %find where there is a big gap between the frames
   missed_rf=find(drf>ifiredbefore*1.75);
   num_missed = length(missed_rf);
    
   if num_missed>1  %pehaps more than one red frame was missed
       warning([titlestr,': ',num2str(num_missed), ' red frames were missed, will try to fix. Check the raw data.']);
   end
   
   %fill the gap, insert the frame in the interval
   while ~isempty(missed_rf)
        for i=1:length(missed_rf)
            ti=rft(missed_rf(i));
            Data(ti+ifiredbefore,redchannel)=1;
        end
        rft=find(Data(:,redchannel)==1);
        ifired=rft(2:end)-rft(1:end-1);
        missed_rf=find(ifired>ifiredbefore*1.75);
   end
    
   % find the beginning and end of the stimulus: 3 consecutive red
   % frames mark the start and end, the stimulus starts with redframe, and
   % the after stimulus gray screen also starts with the red frame
   
   %for each red frame find the number of red frames in a small 
   % neighborhood before and after that frame
   nbwin=550; %small neighborhood: slightly more than three frames
   datalen=length(Data);
   nrft=length(rft);
   nbbefore=zeros(1,nrft); %how many red frames in the neighborhood before
   nbafter=zeros(1,nrft); %how many red frames in the neighborhood after
   for fi=1:length(rft)
       nbbefore(fi)=sum(Data(max(1,rft(fi)-nbwin):rft(fi),2));
       nbafter(fi)=sum(Data(rft(fi):min(datalen,rft(fi)+nbwin),2));
   end

   nbmaxbefore=find(nbbefore==4,1,'first'); %find the first 4 consecutive frames before the current
   nbmaxafter=find(nbafter==3,1,'last')-1; %find the last 3 consecutive frames after the current
   rf1=nbmaxbefore;
   rfe=nbmaxafter;      
   
   %number of red frames during the stimulus presentation
   rfst_before=sum(Data(rft(rf1)+1:rft(rfe),2));
   % remove repeatead red frames during stimulus
   for ri=rf1:rfe-1
       Data(rft(ri)+1:rft(ri)+200,2)=0;
   end
   %how many duplicated frames were removed
   fdiff=rfst_before-sum(Data(rft(rf1)+1:rft(rfe),2));
   %adjust the index of the last stimulus frame
   rfe=rfe-fdiff;
   
   rft=find(Data(:, redchannel));
   ifired=median(rft(rf1+1:rfe)-rft(rf1:rfe-1));
   frame_duration=ifired/5;
   
   %bkgResp=mean(Data(rfe1:rfee,1));%baseline at the end
   bkgResp=mean(Data(rft(rf1):rft(rfe),1)); %baseline as the mean of the recording
%    Data(:,1)=Data(:,1)-bkgResp;
   
%% using lin interpolation find the timing of every frame
   nfr_act=(rfe-rf1+1)*5; %number of red frames of the stimulus
   allframesst=1:nfr_act;% %id of all frames
   rfst=1:5:nfr_act; %id of red frames
   rft_st=rft(rf1:rfe)'; %timing of red frames as recorded
   %liniar interpolation to find timing of all frames
   frametiming=uint32(round(interp1(rfst,rft_st,allframesst,'linear','extrap')));
   
    %% reconstruct frame array 
    [rf,lf,uf,df,nrep]=reconstruct_lines_scanning_rect_with_pause(logfullname);
    %frame_array_lr is the low-res array of the stimulus
    % frame_array  is in the image space. (0,0) in the fly's space is in the
    % bottom rigt corner.
    
    % dirarr is the array of dot's dirrections
    % dirarr()=+/-1 is the movement to the rigth and left. "right" is the
    % increase of the hor coordinate (which in the fly space is reversed, it is movement to
    % the left)
    % dirarr()=+/-2 is the movement upwards and downwards. "up" is the
    % increase of the vertical coordinate 
    % rf, lf,uf,df are the start and end frames of movement to
    % rigth/left/up/down in the image space before transformations to the
    % fly space
    
    %define resolution of the OF picture
    degrees_per_pixel=10;
    screen_span_deg=140; %stimulus spans 140degrees    
    ncol = round(screen_span_deg/degrees_per_pixel);
    scale= ncol/(length(uf)/nrep);
    nrow=round(scale*length(rf)/nrep);
    
        
    %% plot activity with the onset of different frames
    maxval=max(Data(:,1));
    figure,
    plot(Data(:,1),'k'); hold on;
    
    t=frametiming(rf(1,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);    
    p1=plot(t,vals ,'r'); hold on;
    t=frametiming(rf(2,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);    
    p11=plot(t,vals, 'r:'); hold on;
    
    t=frametiming(lf(1,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);  
    p2=plot(t,vals,'b'); hold on;
    t=frametiming(lf(2,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);  
    p22=plot(t,vals,'b:'); hold on;
    
    t=frametiming(uf(1,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);  
    p3=plot(t,vals,'c'); hold on;
    t=frametiming(uf(2,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);  
    p33=plot(t,vals,'c:'); hold on;    
    
    t=frametiming(df(1,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);  
    p4=plot(t,vals,'m'); hold on;
    t=frametiming(df(2,:));
    len=length(t);
    t=repelem(t,3);
    vals=repmat([0, 1*maxval, 0]',len,1);  
    p44=plot(t,vals,'m:'); hold off;
    
    %in the fly space left and rigth, up and down are reversed    
    legend([p1 p2 p3 p4],'Left','Right', 'Up', 'Down'); 
    figname =fullfile(rootDir,[resname,'_data.fig']);
%     savefig(figname);
        
    incr_h_start= frametiming(rf(1,:));
    incr_h_stop = frametiming(rf(2,:));
    decr_h_start= frametiming(lf(1,:));
    decr_h_stop = frametiming(lf(2,:));
    
    incr_v_start = frametiming(uf(1,:));
    incr_v_stop = frametiming(uf(2,:));
    decr_v_start = frametiming(df(1,:));
    decr_v_stop = frametiming(df(2,:));
    
    hor_duration = incr_h_stop(1)-incr_h_start(1);
    hor_lines=size(incr_h_start,2);
    ver_duration = incr_v_stop(1)-incr_v_start(1);
    ver_lines=size(incr_v_start,2);  
    
%     %% uncomment this part to remove subtraction of local baseline    
%     OF_h_incr = cell2mat(arrayfun(@(i) Data(incr_h_start(i):incr_h_start(i)+hor_duration,1)',1:hor_lines, 'UniformOutput', false)');
%     OF_h_decr = cell2mat(arrayfun(@(i) Data(decr_h_start(i):decr_h_start(i)+hor_duration,1)',1:hor_lines, 'UniformOutput', false)');
%     
%     OF_v_incr = cell2mat(arrayfun(@(i) Data(incr_v_start(i):incr_v_start(i)+ver_duration,1),1:ver_lines, 'UniformOutput', false));
%     OF_v_decr = cell2mat(arrayfun(@(i) Data(decr_v_start(i):decr_v_start(i)+ver_duration,1),1:ver_lines, 'UniformOutput', false));
    
    %% compute local baseline activity during "freeze" frames
    % concatinate all lines and sort by the starting time
    alllines = [rf'; lf'; uf'; df'];
    alllines = sortrows(alllines);
    pauses=alllines(2:end,1)-alllines(1:end-1,2);
    freezes_inds=find(pauses>1);
    startfreeze=alllines(freezes_inds,2)+1;
    start_freeze_timing = frametiming(startfreeze);
    freeze_duration=(pauses(freezes_inds(1))-1)*frame_duration;
    
    hor_lines_unique = hor_lines/nrep;
    ver_lines_unique = ver_lines/nrep;
    OF_h_incr=zeros(hor_lines_unique,hor_duration+1);
    OF_h_decr=zeros(hor_lines_unique,hor_duration+1);
    OF_v_incr=zeros(ver_duration+1,ver_lines_unique);
    OF_v_decr=zeros(ver_duration+1,ver_lines_unique);
    
    %change this parameter if you know that the neuron has delay in
    %activity
    response_offset = 0; %0.05*fs; %50ms of the response offset
    %correlations
    for i=1:hor_lines
        %local baseline
        % find the start of the closest freeze
        [mind, pauseind] = min(start_freeze_timing- incr_h_start(i));
        pauset=start_freeze_timing(pauseind);
        bl = mean(Data(pauset:pauset+freeze_duration,1));
        
        ii=mod(i,hor_lines_unique);
        if ii==0; ii=hor_lines_unique; end
        OF_h_incr(ii,:) = OF_h_incr(ii,:)+(Data(incr_h_start(i)+response_offset:incr_h_start(i)+response_offset+hor_duration,1)-bl)';
        OF_h_decr(ii,:) = OF_h_decr(ii,:)+(Data(decr_h_start(i)+response_offset:decr_h_start(i)+response_offset+hor_duration,1)-bl)';        
    end
    
    for i=1:ver_lines
        %local baseline
        % find the start of the closest freeze
        [mind, pauseind] = min(start_freeze_timing- incr_v_start(i));
        pauset=start_freeze_timing(pauseind);
        bl = mean(Data(pauset:pauset+freeze_duration,1));
        
        ii=mod(i,ver_lines_unique);
        if ii==0; ii=ver_lines_unique; end
        OF_v_incr(:,ii) = OF_v_incr(:,ii)+Data(incr_v_start(i)+response_offset:incr_v_start(i)+response_offset+ver_duration,1)-bl;
        OF_v_decr(:,ii) = OF_v_decr(:,ii)+Data(decr_v_start(i)+response_offset:decr_v_start(i)+response_offset+ver_duration,1)-bl;
    end
    OF_h_incr=OF_h_incr./nrep;
    OF_h_decr=OF_h_decr./nrep;
    OF_v_incr=OF_v_incr./nrep;
    OF_v_decr=OF_v_decr./nrep;
    
    %change orientation of the decrease direction
    OF_h_decr = fliplr(OF_h_decr);
    OF_v_decr = flipud(OF_v_decr);
    
    %resize directional RF
    OF_h_incr=imresize(OF_h_incr,[nrow,ncol]);
    OF_h_decr=imresize(OF_h_decr,[nrow,ncol]);
    OF_v_incr=imresize(OF_v_incr,[nrow,ncol]);
    OF_v_decr=imresize(OF_v_decr,[nrow,ncol]);    
    
    OF_hor = OF_h_incr-OF_h_decr;
    OF_ver = OF_v_incr-OF_v_decr;
    %mean RF
    RF= (OF_h_incr + OF_h_decr + OF_v_incr + OF_v_decr)./4;   
    
    %% visualize the results: RF and OF 
    
    maxResp_hi=max(OF_h_incr(:));
    maxResp_hd=max(OF_h_decr(:));
    maxResp_vi=max(OF_v_incr(:));
    maxResp_vd=max(OF_v_decr(:));
    maxResp=max([maxResp_hi, maxResp_hd, maxResp_vi, maxResp_vd]);
    
    minResp_hi=min(OF_h_incr(:));
    minResp_hd=min(OF_h_decr(:));
    minResp_vi=min(OF_v_incr(:));
    minResp_vd=min(OF_v_decr(:));
    minResp=min([minResp_hi, minResp_hd, minResp_vi, minResp_vd]);      
   
    %create one figure with to plots: RF and OF
    fig=figure('Position',[100,100,600,300]);  
    subplot_tight(1,2,1,[0.005,0.05]); 
    lengths=sqrt(OF_hor.^2.+OF_ver.^2);
    lenmax=max(lengths(:));
    OF_hor=OF_hor./lenmax;
    OF_ver=OF_ver./lenmax;
    [x,y] = meshgrid(1:ncol, 1:nrow);
    quiver(x,y,OF_hor,OF_ver, 0,'k','LineWidth',1);     
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);  
    
    subplot_tight(1,2,2,[0.005,0.05]);  
    minRF=min(RF(:));
    maxRF=max(RF(:));
    imagesc(RF,[minRF, maxRF]); 
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    axis equal;
    axis tight;  
    colorbar('northoutside');
    
    ht = suptitle(['OF and RF for ',titlestr]);   
    set(ht,'fontsize',12);
    RF_name=fullfile(rootDir,[resname,'_OF_RF.png']);
%     saveas(fig,RF_name);    
    
    %% visualize the results: RF and OF of each component separately   
      
    fig=figure('Position',[100,100,1200,500]); 
    maxabsResp=max(abs(minResp),abs(maxResp));
    hones=zeros(size(x));
    subplot_tight(2,4,1,[0.005,0.05]);  
    quiver(x,y,OF_h_incr./maxabsResp,hones, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    title('H increase')
    
    subplot_tight(2,4,2,[0.005,0.05]);   
    quiver(x,y,-OF_h_decr./maxabsResp,hones, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    title('H decrease');
    
    subplot_tight(2,4,3,[0.005,0.05]);   
    quiver(x,y, hones, OF_v_incr./maxabsResp, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    title('V increase');
    
    subplot_tight(2,4,4,[0.005,0.05]);   
    quiver(x,y, hones, -OF_v_decr./maxabsResp, 0, 'k','LineWidth',1);       
    set(gca,'XDir','reverse'); 
    axis equal;
    xlim([0 ncol+1]);
    ylim([0 nrow+1]);
    title('V decrease');    
    
    subplot_tight(2,4,5,[0.005,0.05]); 
    imagesc(OF_h_incr,[minResp, maxResp]); 
    set(gca,'YDir','normal'); %(0,0) in fly view is at the bottom right
    set(gca,'Xdir','reverse');
    axis equal;
    axis tight;   
    
    subplot_tight(2,4,6,[0.005,0.05]);    
    imagesc(OF_h_decr,[minResp, maxResp]); 
    set(gca,'YDir','normal');  
    set(gca,'Xdir','reverse');
    axis equal;
    axis tight;
    
    subplot_tight(2,4,7,[0.005,0.05]); 
    imagesc(OF_v_incr,[minResp, maxResp]); %need to flip left right for the fly's space
    set(gca,'YDir','normal');    
    set(gca,'Xdir','reverse');
    axis equal;
    axis tight;
    
    subplot_tight(2,4,8,[0.005,0.05]);       
    imagesc(OF_v_decr,[minResp, maxResp]); %need to flip left right for the fly's space
    set(gca,'Xdir','reverse');
    set(gca,'YDir','normal');
    pos1=get(gca,'Position');
    axis equal;
    axis tight;
   
    cbar_handle =  colorbar('westoutside');   
    set(cbar_handle, 'YAxisLocation','right'); 
    pos21=pos1(1)+pos1(3)+0.01;
    pos22=0.125;%pos1(2)*10;
    pos23=0.01;
    pos24=pos1(4)/2;
    set(cbar_handle,'position',[pos21, pos22,pos23,pos24]);
    suptitle(['Individual components of OF and RF for ',titlestr]);       
    RF_name=fullfile(rootDir,[resname,'_ind_components_RFOF.png']);
%     saveas(fig,RF_name);      
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

