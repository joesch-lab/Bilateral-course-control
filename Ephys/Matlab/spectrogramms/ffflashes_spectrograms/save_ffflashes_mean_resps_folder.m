function save_ffflashes_mean_resps_folder(curfolder)
% for all full field flashes folder get the part of the traces
% correspondign to on and off flash
Fs = 10000; %sampling rate
flash_dur_s=2;
flash_resp=0.5;
tbefore_s=1;
tafter_s=1;

trace_dur=2*round(flash_resp*Fs); %0.5sec before and during the flash
flash_onset = round(flash_resp*Fs)+1;
%the window during which the mean response will be computed
resp_window=round(0.5*Fs);


% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\VS\VS1_4\210128';
    %% get strain, cell type and date from the folder path
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(curfolder);
    
    %% get all stimulus files
    filelist = dir(curfolder);%get list of files
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    stim_name='_full_field_flashes_';
    
    %go thru the list and collect files with'_full_field_flashes_' 
    stim_files={};
    i=1;
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind)
            stim_files{i}=filelist(li).name;            
            i=i+1;
        end
    end
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    total_nrep=0;
    total_nrep_off=0;
    %% iterate thru stim files, get raw ON/OFF traces
    nfiles=length(stim_files);
    ncount=0;
    logiles={};
    prfiles={};
    for i=1:nfiles        
        fulllog_file=fullfile(curfolder,stim_files{i});
                
        rawdata_i=trace_on_of_flashes_one_withpreflashtrace_recording(fulllog_file, flash_dur_s, tbefore_s, tafter_s, flash_resp); 
        if isempty(rawdata_i)            
            continue;
        end
        ncount=ncount+1;
        %structure to pool all data        
        logiles{ncount}=stim_files{i};        
        prfiles{ncount}=rawdata_i.prname;
        data(ncount).trace=rawdata_i.trace;                
        total_nrep=total_nrep+size(data(ncount).trace.on_traces,1);
        total_nrep_off=total_nrep_off+size(data(ncount).trace.off_traces,1);
    end

    if isempty(logiles)
        return;
    end
    
    %% mean and std flash responses
    trace_on=zeros(total_nrep,trace_dur);    
    trace_off=zeros(total_nrep_off,trace_dur);    
    
    %contribution of parts of spectrum: low, mid, total
    resp_on = zeros(total_nrep,1);
    resp_off = zeros(total_nrep_off,1);
        
    %% for each rep and recording make a powerspectrum and average
    nall_on=1;
    nall_off=1;
    for i=1:ncount
        nrepi = size(data(i).trace.on_traces,1);
        for re=1:nrepi 
            trace_on(nall_on,:)=data(i).trace.on_traces(re,:);
            resp_end=min(trace_dur,flash_onset+resp_window);
            %mean response within the window minus the baseline before the
            %flash onset
            %take average of the whole stimulus period and subtract from
            %the prestimulsu period
            resp_on(nall_on)=mean(trace_on(nall_on,flash_onset:resp_end))-mean(trace_on(nall_on,1:flash_onset-1));
            nall_on=nall_on+1;
        end

        nrepi = size(data(i).trace.off_traces,1);
        for re=1:nrepi   
            trace_off(nall_off,:)=data(i).trace.off_traces(re,:);
            resp_end=min(trace_dur,flash_onset+resp_window);
            %mean response within the window minus the baseline before the
            %flash onset
            resp_off(nall_off)=mean(trace_off(nall_off,flash_onset:resp_end))-mean(trace_off(nall_off,1:flash_onset-1));
            nall_off=nall_off+1;
        end 
    end

    %save averages and std    
    exp_data.trace.trace_on_av=mean(trace_on,1);
    exp_data.trace.trace_off_av=mean(trace_off,1);
    exp_data.trace.trace_on_std=std(trace_on,[],1);
    exp_data.trace.trace_off_std=std(trace_off,[],1);
    exp_data.trace.trace_on_all=trace_on;
    exp_data.trace.trace_off_all=trace_off;

    exp_data.resp.resp_on_av=mean(resp_on,1);
    exp_data.resp.resp_off_av=mean(resp_off,1);
    exp_data.resp.resp_on_std=std(resp_on,[],1);
    exp_data.resp.resp_off_std=std(resp_off,[],1);
    
    exp_data.nrep_on= total_nrep;
    exp_data.nrep_off= total_nrep_off;
    exp_data.flash_onset=flash_onset;
    exp_data.flash_chunk=flash_resp;
        
    %create data structure    
    exp_data.folder=curfolder;
    exp_data.logfiles=logiles;
    exp_data.prfiles=prfiles;    
    %cell info
    exp_data.cellinfo.strain = strain;
    exp_data.cellinfo.cellgroup = cellgroup;
    exp_data.cellinfo.celltype = celltype;
    exp_data.cellinfo.datestr = datestr;
    exp_data.cellinfo.cellid = cellid;       

    %%to check the data plot raw extracted traces during and off flash
%     figure; t=tiledlayout(2,1);
%     xvals=1:trace_dur;
%     xvals_plot=xvals/Fs;
%     flash_onset_plot=flash_onset/Fs;
% 
%     nexttile;
%     plot(xvals_plot, trace_on',"Color",[0.5,0.5,0.5]);
%     hold on;
%     plot(xvals_plot, exp_data.trace.trace_on_av,"Color",[0,0,0],"LineWidth",1);
%     hold on;
%     maxval=max(exp_data.trace.trace_on_av);
%     minval=min(exp_data.trace.trace_on_av);
%     plot([flash_onset_plot,flash_onset_plot],[minval,maxval],'b');
%     ylabel('ON');
% 
%     nexttile;
%     plot(xvals_plot, trace_off',"Color",[0.5,0.5,0.5]);
%     hold on;
%     plot(xvals_plot, exp_data.trace.trace_off_av,"Color",[0,0,0],"LineWidth",1);
%     hold on;
%     maxval=max(exp_data.trace.trace_off_av);
%     minval=min(exp_data.trace.trace_off_av);
%     plot([flash_onset_plot,flash_onset_plot],[minval,maxval],'b');
%     ylabel('OFF');
%     title(t,curfolder,'Interpreter','none');
    
    %save    
    if isempty(cellid)
        filename_i=['flashes_on_off_traces_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['flashes_on_off_traces_',strain,'_',celltype,'_',datestr,'_',cellid];        
    end    
    fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
    save(fullfile_i,'exp_data','-v7.3');  
end

