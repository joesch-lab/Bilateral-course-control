% given a folder, load all stimuli_grating_stimuli_XXXXXX files with 
% grating in different directions
% for each set of params (speed, spatial_frequency) average recordings
% each set of params should be stored in different files

function save_grating_stimuli_raw_data_folder(curfolder)  
    allgratings=['*stimuli_grating_stimuli_*.log'];
    stim_name='grating_stimuli_[\d{4}]';  

    filelist = dir(fullfile(curfolder, allgratings));%get list of grating files
    
    %go thru the list and collect all folders with'grating_stimuli_' 
    files2analyze={};
    i=1;
    for li =1:length(filelist)    
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) 
            files2analyze{i}=filelist(li).name;            
            i=i+1;  
        end
    end

    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    nfiles=length(files2analyze);
    set_params=[];
    for i=1:nfiles
        logfullname=fullfile(curfolder,files2analyze{i});
        stimparam = load_grating_stimparam(logfullname);
        set_params=[set_params;[stimparam.speed, stimparam.sp_freq,{files2analyze{i}}]];
    end

    set_params_sub=cell2mat(set_params(:,1:2));
    unique_paramset=unique(set_params_sub,'rows');
    nsets=size(unique_paramset,1);

    %for a set of params load the files and average
    for si=1:nsets
        %find files with these params        
        fids=find(ismember(set_params_sub,unique_paramset(si,:),"rows"));
        
        exp_data={};
        exp_data.stimparam.speed=unique_paramset(si,1);
        exp_data.stimparam.sp_freq=unique_paramset(si,2);
        exp_data.logfiles={};
        exp_data.prfiles={};

        for fi=1:length(fids)
            logfile=set_params{fids(fi),3};
            logfullname=fullfile(curfolder,logfile);
            one_rec_data = grating_stimuli_one_recording(logfullname);
            if isempty(one_rec_data)
                continue;
            end
            exp_data.logfiles{fi}=logfile;
            exp_data.prfiles{fi}=one_rec_data.prname;
    
            if isfield(exp_data, 'traces')==0
                exp_data.cellinfo=one_rec_data.cellinfo;

                exp_data.traces=one_rec_data.all_repsp;
                exp_data.baseline=one_rec_data.all_baseline;
                exp_data.mean_baseline = one_rec_data.mean_baseline; 
                exp_data.stimparam.dir_vec=0:45:359;
            else
                rec_len=min(size(exp_data.traces,3),size(one_rec_data.all_repsp,3));
                rec_len_bl=min(size(exp_data.baseline,3),size(one_rec_data.all_baseline,3));
                %trim traces to min length
                if size(exp_data.traces,3)>rec_len
                    exp_data.traces(:,:,rec_len+1:end)=[];
                end
                if size(one_rec_data.all_repsp,3)>rec_len
                    one_rec_data.all_repsp(:,:,rec_len+1:end)=[];
                end
                %trim baseline to min length
                if size(exp_data.baseline,3)>rec_len_bl
                    exp_data.baseline(:,:,rec_len_bl+1:end)=[];
                end
                if size(one_rec_data.all_baseline,3)>rec_len_bl
                    one_rec_data.all_baseline(:,:,rec_len_bl+1:end)=[];
                end
                exp_data.traces=cat(2,exp_data.traces,one_rec_data.all_repsp);
                exp_data.baseline=cat(2,exp_data.baseline,one_rec_data.all_baseline);
                exp_data.mean_baseline = cat(2,exp_data.mean_baseline,one_rec_data.mean_baseline); 
                
            end
        end
        %% extract spikes and add to the file
        rawspikes=zeros(size(exp_data.traces));
        baseline_spikes=zeros(size(exp_data.baseline));
        for j=1:size(exp_data.traces,2)
            [raw_spike_data, raw_all_spikes]=get_spike_data(squeeze(exp_data.traces(:, j, :)));
            rawspikes(:, j, :)=raw_spike_data;
            [baseline_spike_data, baseline_all_spikes]=get_spike_data(squeeze(exp_data.baseline(:, j, :)));
            baseline_spikes(:, j, :)=baseline_spike_data;
        end
        %%
        spikerate_rep=sum(rawspikes, 3)*(10000/size(rawspikes, 3));
        baseline_spikerate_rep=sum(baseline_spikes, 3)*(10000/size(baseline_spikes, 3));
        exp_data.meanspike_rate=mean(spikerate_rep - baseline_spikerate_rep, 2);
        exp_data.baselinespike_rate=mean(baseline_spikerate_rep, 2); 

        %average traces
        exp_data.traces_av=squeeze(mean(exp_data.traces,2,'omitnan'));
        exp_data.baseline_av=squeeze(mean(exp_data.baseline,2,'omitnan'));
        exp_data.traces_std=squeeze(std(exp_data.traces,[],2,'omitnan'));
        exp_data.baseline_std=squeeze(std(exp_data.baseline,[],2,'omitnan'));
        exp_data.nrep=size(exp_data.traces,2);

        %save data
        
        spfreq=num2str(round(exp_data.stimparam.sp_freq,2));
        tempstr=strsplit(spfreq,'.');
        spfreq=['spf',tempstr{2}];
        speed_str=num2str(round(exp_data.stimparam.speed));

        if isempty(exp_data.cellinfo.cellid)
            filename_i=['grating_folder_',exp_data.cellinfo.strain_type,'_',exp_data.cellinfo.celltype,'_',exp_data.cellinfo.datestr,'_dx',speed_str,'_',spfreq];        
        else
            filename_i=['grating_folder_',exp_data.cellinfo.strain_type,'_',exp_data.cellinfo.celltype,'_',exp_data.cellinfo.datestr,'_',exp_data.cellinfo.cellid,'dx',speed_str,'_',spfreq];        
        end
        fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
        save(fullfile_i,'exp_data','-v7.3');
    end
end


function stimparam = load_grating_stimparam(logfullname)
% read the file with stimuli parameters
    f=fopen(logfullname);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'grating_stimuli'))
        disp('The log file does not contain information about grating stimuli.');
        fclose(f);
        return;
    end

    tline=fgetl(f); %start time of the recording    

    tline=fgetl(f); %end time of the recording
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparam.time=datevec(tstr,formatIn);
    
    tline=fgetl(f); %time per grating 
    tline=fgetl(f); %Flip rate: 
    
    tline=fgetl(f);
    tstr=tline(length('Repetitions: ')+1:end);    
    stimparam.nrep=str2double(tstr); 
    
    tline=fgetl(f);
    tstr=tline(length('Speed: ')+1:end);    
    stimparam.speed=str2double(tstr); 
        
    tline=fgetl(f);
    tstr=tline(length('Baseline recording: ')+1:end);     
    stimparam.baseline_dur=uint32(str2num(tstr)); 
    
    tline=fgetl(f);
    tstr=tline(length('Spatial frequency: ')+1:end);     
    stimparam.sp_freq=str2double(tstr);     
    
    fclose(f);
end
