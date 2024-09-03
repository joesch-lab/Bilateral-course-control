function save_velocity_tuning_raw_data_folder(curfolder)
% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\VS\VS1_4\210128';
    %% get all stimulus files
    filelist = dir(curfolder);%get list of files
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    stim_name='S_Vel_';
    rec_duration=4*10000;
    %go thru the list and collect files with'grating_stimuli_*S_Vel_' 
    stim_files={};
    i=1;
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'.log','IgnoreCase', true)  
            stim_files{i}=filelist(li).name;            
            i=i+1;
        end
    end
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    %% get all possible speeds for data grouping
    [sp_arr,~] = all_speed_spfreq_param_vel_tuning(curfolder);   
    
    %% iterate thru stim files, group files with the same param, pool data
    nfiles=length(stim_files);
    file_analized=zeros(1,nfiles);
    for i=1:nfiles
        if file_analized(i)==1
            continue;
        end
        ncount=1;
        fulllog_file=fullfile(curfolder,stim_files{i});
        stimparam_i=read_velocity_tuning_stimuli_log(fulllog_file);
        %get all combinations of speed and directions
        if stimparam_i.isHS==1
            dir_arr=[0,180];
        else
            dir_arr=[90,270];
        end
         %get all combinations of speed and directions
        [x,y]=meshgrid(sp_arr, dir_arr);
        x=x(:); y=y(:);
        speed_dir=[x';y']';
        nspeedsdir=size(speed_dir,1); 
        rawtraces=NaN(nspeedsdir,1,rec_duration);
        baseline_rawtraces=NaN(nspeedsdir,1,5000);

        rawdata_i=vel_tuning_one_recording(fulllog_file); 
        if isempty(rawdata_i)
            file_analized(i)=1;
            continue;
        end
        %structure to pool all data
        logiles={};
        prfiles={};
        logiles{ncount}=stim_files{i};        
        prfiles{ncount}=rawdata_i.prname;
        nreps=rawdata_i.nrep;
        nrows=size(rawdata_i.raw_data,2);        
        trace=reshape([rawdata_i.raw_data(:).trace],[],nrows)';
        baseline_trace=reshape([rawdata_i.raw_data(:).baseline_trace],[],nrows)';
        reclen=length(trace);
        minlen=reclen;
        %corresponding indices of speed entries
        speedsi = [rawdata_i.raw_data.speed];
        diri=[rawdata_i.raw_data.dir];
        speed_dir_i=[speedsi;diri]';
        
        rawspikes=NaN(nspeedsdir,nreps,minlen);
        baseline_spikes=NaN(nspeedsdir,nreps,5000);

        inds = zeros(size(speed_dir_i,1),1);
        for si=1:length(inds)
            indsi = find(ismember(speed_dir, speed_dir_i(si,:),'rows'));
            inds(si)=indsi;
        end    
    
        %find indices in the sorted speed array        
        rawtraces(inds,1:rawdata_i.nrep,1:reclen)=trace;
        baseline_rawtraces(inds,1:rawdata_i.nrep,:)=baseline_trace(:, 1:5000);
        %add spike data
        %
        file_analized(i)=1;
        %collect all recordings with the same params
        for j=i+1:nfiles
            if file_analized(j)==1
                continue;
            end
            fulllog_file_j=fullfile(curfolder,stim_files{j});
            stimparam_j=read_velocity_tuning_stimuli_log(fulllog_file_j);
            %check params that matter for averaging:
            % cell type and spatial frequency
            if stimparam_i.isHS~=stimparam_j.isHS || stimparam_i.sp_freq~=stimparam_j.sp_freq
                continue;
            end
            rawdata_j=vel_tuning_one_recording(fulllog_file_j); 
            if isempty(rawdata_j)
                file_analized(j)=1;
                continue;
            end
            %structure to pool all data
            logiles{end+1}=stim_files{j};      
            prfiles{end+1}=rawdata_j.prname;
            nreps=nreps+rawdata_j.nrep;
            
            lastrepind=size(rawtraces,2);
            
            nrows=size(rawdata_j.raw_data,2);        
            trace=reshape([rawdata_j.raw_data(:).trace],[],nrows)';
            baseline_trace=reshape([rawdata_i.raw_data(:).baseline_trace],[],nrows)';
            reclen=length(trace);
            minlen=min(minlen,reclen);
            %corresponding indices of speed entries
            speedsi = [rawdata_j.raw_data.speed];
            diri=[rawdata_j.raw_data.dir];
            speed_dir_i=[speedsi;diri]';

            inds = zeros(size(speed_dir_i,1),1);
            for si=1:length(inds)
                indsi = find(ismember(speed_dir, speed_dir_i(si,:),'rows'));
                inds(si)=indsi;
            end               
            rawtraces(inds,lastrepind+1:lastrepind+rawdata_j.nrep,1:reclen)=trace;
            baseline_rawtraces(inds,lastrepind+1:lastrepind+rawdata_j.nrep,:)=baseline_trace(:, 1:5000);
            file_analized(j)=1;            
        end                
        
        %crop the time to the smallest recording
        rawtraces=rawtraces(:,:,1:minlen);
        %%
        for j=1:nreps
            [raw_spike_data, raw_all_spikes]=get_spike_data(squeeze(rawtraces(:, j, :)));
            rawspikes(:, j, :)=raw_spike_data;
            [baseline_spike_data, baseline_all_spikes]=get_spike_data(squeeze(baseline_rawtraces(:, j, :)));
            baseline_spikes(:, j, :)=baseline_spike_data;
        end
        %%
        %average reps
        meanval_rep=mean(mean(rawtraces,3),2);
        std_rep=std(mean(rawtraces,3),[],2);
        spikerate_rep=sum(rawspikes, 3)*(10000/size(rawspikes, 3));
        baseline_spikerate_rep=sum(baseline_spikes, 3)*(10000/size(baseline_spikes, 3));
        meanspike_rate=mean(spikerate_rep - baseline_spikerate_rep, 2);
        baselinespike_rate=mean(baseline_spikerate_rep, 2);
        
        %create data structure
        clear exp_data;
        exp_data.folder=curfolder;
        exp_data.logfiles=logiles;
        exp_data.prfiles=prfiles;
        exp_data.isHS = stimparam_i.isHS;
        exp_data.sp_freq = stimparam_i.sp_freq;
        exp_data.tmove= stimparam_i.tmove;
        exp_data.nrep= nreps;
        exp_data.raw_traces = rawtraces;
        exp_data.baseline_traces = baseline_rawtraces;
        exp_data.mean_vals= meanval_rep;
        exp_data.std_vals= std_rep;
        exp_data.speed_dir= speed_dir;
        exp_data.mean_spikes= meanspike_rate;
        exp_data.baseline_spikes= baselinespike_rate;
        
        %save
        %get strain, cell type and date from the folder path
        [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(curfolder);
        if exp_data.isHS == 1;
            stimHS='stimHS';
        else
            stimHS='stimVS';
        end
        spfreq=num2str(round(exp_data.sp_freq,2));
        tempstr=strsplit(spfreq,'.');
        spfreq=['spf',tempstr{2}];
        if isempty(cellid)
            filename_i=['vel_tuning_',strain,'_',celltype,'_',datestr,'_',stimHS,'_',spfreq];        
        else
            filename_i=['vel_tuning_',strain,'_',celltype,'_',datestr,'_',cellid,'_',stimHS,'_',spfreq];        
        end
        fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
        
%         %check how many files with this name exist
%         nsame_files = count_files_namebase(resfolder,filename_i);
%         if nsame_files==0
%             fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
%         else
%             fullfile_i=fullfile(resfolder,[filename_i,'_',num2str(nsame_files),'.mat']);
%         end
%             
        
        %issue a warning if the file with this name already exist
%         if nsame_files>0
%             warning(['File ',fullfile_i,' already exists. Check why.']);            
%         end
        save(fullfile_i,'exp_data','-v7.3');
    end
  
end


function  nfiles = count_files_namebase(curfolder,filebase)
    filelist = dir(curfolder);%get list of files
    filelist = filelist(~[filelist.isdir]);  %remove folders from list
    
    %go thru the list and count the files containing filebase in the name    
    nfiles=0;
    for li =1:length(filelist)        
        if contains(filelist(li).name,filebase,'IgnoreCase', true)                       
           nfiles=nfiles+1;
        end
    end
end