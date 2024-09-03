function save_contrast_levels_raw_data_folder(curfolder)
% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND DB331\VS\VS1_4\220209\male7\';
    %% get all stimulus files    
    allsubfolders=['**',filesep,'*_stimuli_contrast_levels_*log.mat'];
    stim_files = dir(fullfile(curfolder, allsubfolders));

%     stim_name='_stimuli_contrast_levels_';
%     rec_duration=4*10000;
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    
    %% iterate thru stim files, group files with the same param, pool data
    nfiles=length(stim_files);
    file_analized=zeros(1,nfiles);
    for i=1:nfiles
        if file_analized(i)==1
            continue;
        end
        ncount=1;
        fulllog_file=fullfile(curfolder,stim_files(i).name);
        stimparam_i=load(fulllog_file);
        n_contrast_i = size(stimparam_i.contrast_levels,1);
        if n_contrast_i~=5
            file_analized(i)=1;
            continue;
        end
        dir_i= unique(stimparam_i.allconfigs(:,2));
        speed_i=stimparam_i.dx;
        sp_f_i = stimparam_i.freqCyclesPerPix;
        nconfigs = size(stimparam_i.allconfigs,1);
        nrep = stimparam_i.repetitions;

        rawdata_i=contrast_levels_one_recording(fulllog_file); 
        rec_duration=size(rawdata_i.trace,2);

        %get all combinations of speed and directions
        rawtraces=NaN(nconfigs,nrep,rec_duration);
        baselines=NaN(nconfigs,nrep);
        baseline_traces=NaN(nconfigs,nrep,5000);
        rawspikes=NaN(nconfigs,nrep,rec_duration);
        baseline_spikes=NaN(nconfigs,nrep,5000);

        if isempty(rawdata_i)
            file_analized(i)=1;
            continue;
        end
        %structure to pool all data
        logiles={};
        prfiles={};
        logiles{ncount}=stim_files(i).name;        
        prfiles{ncount}=rawdata_i.prname;
        nreps=nrep;               
       
        minlen=rec_duration;
        %corresponding indices of contrast and direction entries        
        inds = stimparam_i.allconfigs_reps;
        
        
        %add reps
        for ri=1:nrep
            ind1=(ri-1)*nconfigs+1;
            indn=ri*nconfigs;            
            disp(inds(ind1:indn));
            rawtraces(inds(ind1:indn),ri,1:rec_duration)=rawdata_i.trace(ind1:indn,:);
            baselines(inds(ind1:indn),ri)=rawdata_i.baseline(ind1:indn);
            baseline_traces(inds(ind1:indn),ri,1:5000)=rawdata_i.baseline_traces(ind1:indn,:);
        end

        file_analized(i)=1;
        %collect all recordings with the same params
        for j=i+1:nfiles
            if file_analized(j)==1
                continue;
            end
            fulllog_file_j=fullfile(curfolder,stim_files(j).name);
            stimparam_j=load(fulllog_file);
            n_contrast_j = size(stimparam_j.contrast_levels,1);
            if n_contrast_j~=5
                file_analized(j)=1;
                continue;
            end
            %check params that matter for averaging:            
            dir_j= unique(stimparam_j.allconfigs(:,2));
            speed_j=stimparam_j.dx;
            sp_f_j = stimparam_j.freqCyclesPerPix;
            if ~all(dir_i==dir_j) || speed_i~=speed_j || sp_f_i~=sp_f_j
                continue;
            end
            
            rawdata_j=contrast_levels_one_recording(fulllog_file_j); 
            if isempty(rawdata_j)
                file_analized(j)=1;
                continue;
            end
            %structure to pool all data
            logiles{end+1}=stim_files(j).name;      
            prfiles{end+1}=rawdata_j.prname;
            nrep=stimparam_j.repetitions;
            nreps=nreps+nrep;
            
            lastrepind=size(rawtraces,2);            
            rec_duration=size(rawdata_j.trace,2);                             
            minlen=min(minlen,rec_duration);
            %corresponding indices of contrast and directions
            inds = stimparam_j.allconfigs_reps;                           
            %add reps
            for ri=1:nrep
                ind1=(ri-1)*nconfigs+1;
                indn=ri*nconfigs;            
                rawtraces(inds(ind1:indn),lastrepind+ri,1:minlen)=rawdata_j.trace(ind1:indn,1:minlen);
                baselines(inds(ind1:indn),lastrepind+ri)=rawdata_j.baseline(ind1:indn);
                baseline_traces(inds(ind1:indn),lastrepind+ri,1:5000)=rawdata_i.baseline_traces(ind1:indn,:);
            end
            file_analized(j)=1;
        end                
        
        %crop the time to the smallest recording
        rawtraces=rawtraces(:,:,1:minlen);
%         baseline_traces=baseline_traces(:,:,end-5000:end);

        %% extract spikes and add to the file
        for j=1:nreps
            [raw_spike_data, raw_all_spikes]=get_spike_data(squeeze(rawtraces(:, j, :)));
            rawspikes(:, j, :)=raw_spike_data;
            [baseline_spike_data, baseline_all_spikes]=get_spike_data(squeeze(baseline_traces(:, j, :)));
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
        exp_data.contrast_levels=stimparam_i.contrast_levels;
        exp_data.allconfigs=stimparam_i.allconfigs;
        exp_data.speed = speed_i;
        exp_data.sp_freq = sp_f_i;
        exp_data.tmove= stimparam_i.t;
        exp_data.nrep= nreps;
        exp_data.raw_traces = rawtraces;
        exp_data.baselines = baselines;
        exp_data.baseline_traces = baseline_traces;
        exp_data.mean_vals= meanval_rep;
        exp_data.std_vals= std_rep;
        exp_data.mean_spikes= meanspike_rate;
        exp_data.baseline_spikes= baselinespike_rate;
        %save
        %get strain, cell type and date from the folder path
        [strain,~,celltype,datestr,cellid] = cell_info_from_path(curfolder);
        
        spfreq=num2str(round(exp_data.sp_freq,2));
        tempstr=strsplit(spfreq,'.');
        spfreq=['spf',tempstr{2}];
        dxstr=['dx',num2str(speed_i)];
        if isempty(cellid)
            filename_i=['grating_contrast_levels_',strain,'_',celltype,'_',datestr,'_',dxstr,'_',spfreq];        
        else
            filename_i=['grating_contrast_levels_',strain,'_',celltype,'_',datestr,'_',cellid,'_',dxstr,'_',spfreq];        
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
%         
%         %issue a warning if the file with this name already exist
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