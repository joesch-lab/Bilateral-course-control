% %the repo has the following hierarchy:
% -- strain
%     -- cell type
%         --date
%             --animal
%             
% if there are no subfolders in the date folder, it means there was only 
% one animal per day
%     The script will iterate thru the subfolders, collect and average 
% recordings per animal. The result will be stored in the folder 
% repo/strain/celltype/grating/date_animal_av.mat

function average_grating_per_animal(repo_folder,celltype, speed, sp_freq)    
    %get all subfolders - all strain types
    filelist = dir(fullfile(repo_folder));
    filelist =  filelist([filelist.isdir]);
    exclude_patterns=[".","res"];
    filelist = exclude_subfolders(filelist, exclude_patterns);
    
    
    %iterate strain types
    n_strains=length(filelist);
    for si =1:n_strains
        curr_strain=filelist(si).name;
        strain_folder=fullfile(filelist(si).folder,curr_strain);
        cellgrp_list = dir(strain_folder);
        cellgrp_list =  cellgrp_list([cellgrp_list.isdir]);
        %exclude results and links to .. folders
        cellgrp_list = exclude_subfolders(cellgrp_list, exclude_patterns);
        n_cellgrptypes=length(cellgrp_list);
        for cgi =1:n_cellgrptypes
            curr_cellgrp=cellgrp_list(cgi).name;
            cellgrp_folder=fullfile(strain_folder,curr_cellgrp);     
            celltype_list = dir(cellgrp_folder);
            celltype_list =  celltype_list([celltype_list.isdir]);
            celltype_list = exclude_subfolders(celltype_list, exclude_patterns);

            %iterate thru cell types 
            n_celltypes=length(celltype_list);
            for ci =1:n_celltypes
                curr_cell=celltype_list(ci).name;
                if ~contains(curr_cell, celltype)
                    continue;
                end                
                cell_folder=fullfile(cellgrp_folder, curr_cell);
                date_list = dir(cell_folder);
                date_list =  date_list([date_list.isdir]);
                %exclude results and links to .. folders
                date_list = exclude_subfolders(date_list, exclude_patterns);
                
                %create grating res folder
                gr_folder=fullfile(cell_folder,'grating_stimuli_res');
                if ~exist(gr_folder,'dir')
                    mkdir(gr_folder)
                end

                n_days=length(date_list);
                %iterate thru dates            
                for di =1:n_days
                    curr_date=date_list(di).name;
                    datefolder=fullfile(cell_folder, curr_date);
                    animal_list = dir(datefolder);
                    animal_list =  animal_list([animal_list.isdir]);
                    %exclude results and links to .. folders
                    animal_list = exclude_subfolders(animal_list, exclude_patterns);

                    n_flies=length(animal_list);
                    %if n_flies==0 there is acutally one fly per that date without
                    %a subfolder structure, otherwise more than one
                    if n_flies==0
                        fly_av_vals = process_grating_stimuli_one_fly_folder(datefolder, speed, sp_freq);
                        fly_av_vals.strain_type=curr_strain;
                        fly_av_vals.cell_type = curr_cell;
                        full_fly_name=[curr_date,'_',fly_av_vals.fly_name];
                        fly_av_vals.fly_name=full_fly_name; 
                        file_ending=['_dx',num2str(speed),'_spf',num2str(sp_freq*100),'.mat'];
                        file_name=fullfile(gr_folder,[full_fly_name,file_ending]);
                        save(file_name,'fly_av_vals');
                    else
                        %iterate per animal
                        for fi =1:n_flies
                            curr_fly=animal_list(di).name;
                            flyfolder = fullfile(datefolder, curr_fly);
                            fly_av_vals = process_grating_stimuli_one_fly_folder(flyfolder, speed, sp_freq);
                            fly_av_vals.strain_type=curr_strain;
                            fly_av_vals.cell_type = curr_cell;
                            full_fly_name=[curr_date,'_',fly_av_vals.fly_name];
                            fly_av_vals.fly_name=full_fly_name; 
                            file_ending=['_dx',num2str(speed),'_spf',num2str(sp_freq*100),'.mat'];
                            file_name=fullfile(gr_folder,[full_fly_name,file_ending]);
                            save(file_name,'fly_av_vals');
                        end            
                    end
                end
            end
        end
    end
end
    
    
%     

function average_vals = process_grating_stimuli_one_fly_folder(fly_folder, speed, sp_freq)
    %for a fly folder average all grating recordings
    %grating stimuli is already processed, the results are stored in the
    %res subfolder
    stim_name='grating_stimuli_[\d{4}]';    
    all_entries={};
    i=1;
    min_baseline_dur=Inf;
    min_resp_dur=Inf;
    res_subfolder = fullfile(fly_folder,'res');
    filelist=dir(res_subfolder);
    filelist=filelist(~[filelist.isdir]);
    nfiles = length(filelist);
    
    for li =1:nfiles        
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'_raw_ds.mat','IgnoreCase', true)              
            matfile=fullfile(filelist(li).folder, filelist(li).name);
            
            %load each file and check the file type
            load(matfile);
            %now the structure exp_data is available            
            if exp_data.speed~=speed
                continue;
            end
            if exp_data.sp_freq~=sp_freq
                continue;
            end
            %for matching files collect activity averaged per recording
            
            %responses to all directions
            all_entries(i).av_resp=squeeze(exp_data.av_responses);
            all_entries(i).av_baseline=squeeze(exp_data.av_baseline);
            resp_dur=size(all_entries(i).av_resp,2);
            baseline_dur=size(all_entries(i).av_baseline,2);
            
            min_baseline_dur=min(min_baseline_dur, baseline_dur);
            min_resp_dur=min(min_resp_dur, resp_dur);  
    
            if i==1  
                animal_name=exp_data.prname;
                splitpos = strfind(animal_name,'_');
                animal_name=animal_name(1:splitpos(end)-1);
                dir_vector=exp_data.dir_vector;
                ndir=length(exp_data.dir_vector);
            end
            i=i+1;
            clear exp_data;
        end
    end
   
    ngr_files=i-1;
    if ngr_files==0
        return;
    end
    %compute average responses per animal
    resp_vals=zeros(ngr_files,ndir,min_resp_dur);
    bl_vals=zeros(ngr_files,ndir,min_baseline_dur);
    for j=1:ngr_files
        resp_vals(j,:,:)=all_entries(j).av_resp(:,1:min_resp_dur);
        bl_vals(j,:,:)=all_entries(j).av_baseline(:,1:min_baseline_dur);
    end
    
    average_vals.raw_resp=resp_vals;
    average_vals.raw_bl=bl_vals;
    average_vals.av_resp=squeeze(mean(resp_vals,1));
    average_vals.std_resp=squeeze(std(resp_vals,[],1));
    average_vals.av_bl=squeeze(mean(bl_vals,1));
    average_vals.std_bl=squeeze(std(bl_vals,[],1));
    average_vals.n_grating_files=ngr_files;   
    average_vals.fly_name=animal_name;
    average_vals.dir_vector=dir_vector;
    average_vals.speed=speed;
    average_vals.sp_freq=sp_freq;
end

function filelist = exclude_subfolders(filelist, exclude_patterns)
    exclude_list=[];
    for li =1:length(filelist)
        if contains(filelist(li).name,exclude_patterns)
            exclude_list=[li,exclude_list];
        end
    end
    filelist(exclude_list)=[];
end
