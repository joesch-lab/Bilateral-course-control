%for a given cell type and stimulus conditions collect all recordings
function collect_cell_type_recordings(parent_folder,celltype, speed, sp_freq)    
     %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*.*'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    %go thru the list, analize each log file with 'grating_stimuli' in the title
    stim_name='grating_stimuli_[\d{4}]';    
    all_entries={};
    i=1;
    min_baseline_dur=Inf;
    min_resp_dur=Inf;
    for li =1:length(filelist)
       
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'_raw_ds.mat','IgnoreCase', true)              
            matfile=fullfile(filelist(li).folder, filelist(li).name);
            celltype_i=get_celltype_from_res_folder_path(matfile);
            if ~strcmp(celltype_i,celltype)
                continue;
            end
            %load each file and check the file type
            load(matfile);
            %now the structure exp_data is available            
            if exp_data.speed~=speed
                continue;
            end
            if exp_data.sp_freq~=sp_freq
                continue;
            end
            %for matching files collect activity averaged per animal
            all_entries(i).strain_type=exp_data.strain_type;
            %responses to all directions
            all_entries(i).av_resp=squeeze(exp_data.av_responses);
            all_entries(i).av_baseline=squeeze(exp_data.av_baseline);
            resp_dur=size(all_entries(i).av_resp,2);
            baseline_dur=size(all_entries(i).av_baseline,2);
            
            min_baseline_dur=min(min_baseline_dur, baseline_dur);
            min_resp_dur=min(min_resp_dur, resp_dur);
  
    
            if i==1
                exp_info.speed=speed;
                exp_info.sp_freq=sp_freq;
                exp_info.cell_type=celltype;
                exp_info.dir_vector=exp_data.dir_vector;
                ndir=length(exp_data.dir_vector);
            end
            i=i+1;
            clear exp_data;
        end
    end
   
    %compute average responses per strain type
    names_table=struct2table(all_entries).strain_type;
    all_strain_types =  unique(names_table);
    nstrains=size(all_strain_types,1);
    average_vals={};
    for i=1:nstrains
        curr_strain=all_strain_types(i);
        average_vals(i).strain_type=curr_strain;
        idk=find(strcmp(names_table, curr_strain)); 
        nflies=length(idk);
        resp_vals=zeros(nflies,ndir,min_resp_dur);
        bl_vals=zeros(nflies,ndir,min_baseline_dur);
        if length(idk)==0
            average_vals(i).nflies=0;
            average_vals(i).av_resp=zeros(ndir,min_resp_dur);
            average_vals(i).std_resp=zeros(ndir,min_resp_dur);
            average_vals(i).av_bl=zeros(ndir,min_baseline_dur);
            average_vals(i).std_bl=zeros(ndir,min_baseline_dur);
            average_vals(i).nflies=nflies;
            continue;
        end
        for j=1:length(idk)
            id=idk(j);
            resp_vals(j,:,:)=all_entries(id).av_resp(:,1:min_resp_dur);
            bl_vals(j,:,:)=all_entries(id).av_baseline(:,1:min_baseline_dur);
        end
        average_vals(i).raw_resp=resp_vals;
        average_vals(i).raw_bl=bl_vals;
        average_vals(i).av_resp=squeeze(mean(resp_vals,1));
        average_vals(i).std_resp=squeeze(std(resp_vals,[],1));
        average_vals(i).av_bl=squeeze(mean(bl_vals,1));
        average_vals(i).std_bl=squeeze(std(bl_vals,[],1));
        average_vals(i).nflies=nflies;
    end
    
    mat_file = [celltype,'_all_recordings.mat'];
    mat_file =fullfile(parent_folder,mat_file);
    save(mat_file,'average_vals', 'exp_info');    
end

function celltype = get_celltype_from_res_folder_path(logfile)
    [filepath,filename,~]=fileparts(logfile);
    [filepath,resfolder,~]=fileparts(filepath);
    [filepath,datestr,~]=fileparts(filepath);
    [filepath,celltype,~]=fileparts(filepath);
end