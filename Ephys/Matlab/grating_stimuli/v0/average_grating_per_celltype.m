% %the repo has the following hierarchy:
% -- strain
%     -- cell type
%         --date
%             --animal
%          --grating_stimuli_res
%               --animal_average_stimuli_results
%             
%     The script will iterate thru all strain folders and for the given 
% celltype will average all recordings the result will be stored in
% repo/grating_stimuli_res/grating_stimuli_strain_celltype.mat

function average_grating_per_celltype(repo_folder,celltype, speed, sp_freq) 
    resfolder = fullfile(repo_folder,'grating_stimuli_res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder)
    end
    
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
                %check if there is a grating_stimuli_res folder
                grresfolder = fullfile(cell_folder,'grating_stimuli_res');
                if ~exist(grresfolder,'dir')
                    continue;
                end                
                %average grating stimuli results accross same cell type
                celltype_av_vals = average_grating_stimuli_celltype(grresfolder,speed,sp_freq);
                if isempty(celltype_av_vals)
                    continue;
                end
                celltype_av_vals.strain_type=curr_strain;
                celltype_av_vals.cell_type = curr_cell;
                celltype_file=['grating_stimuli_', curr_strain,'_',curr_cell];
                file_ending=['_dx',num2str(speed),'_spf',num2str(sp_freq*100),'.mat'];
                file_name=fullfile(resfolder,[celltype_file,file_ending]);
                save(file_name,'celltype_av_vals');   
            end
        end
    end
end
    
    
function average_vals = average_grating_stimuli_celltype(gr_res_files,speed,sp_freq)
    %for the folder with grating results per fly collect all the data and
    %make an average per celltype
    
    all_entries={};
    i=1;
    min_baseline_dur=Inf;
    min_resp_dur=Inf;
    
    filelist=dir(gr_res_files);
    filelist=filelist(~[filelist.isdir]);
    nfiles = length(filelist);
    
    for li =1:nfiles   
        %exclude any file that doesn't start with the date
        if ~startsWith(filelist(li).name,digitsPattern(6))
            continue;
        end
        matfullfile= fullfile(gr_res_files,filelist(li).name);
        load(matfullfile);
        
        %now the structure fly_av_vals is available            
        if fly_av_vals.speed~=speed
            continue;
        end
        if fly_av_vals.sp_freq~=sp_freq
            continue;
        end
        %for matching files collect activity averaged per fly
            
        %responses to all directions
        all_entries(i).av_resp=squeeze(fly_av_vals.av_resp);
        all_entries(i).av_baseline=squeeze(fly_av_vals.av_bl);
        resp_dur=size(all_entries(i).av_resp,2);
        baseline_dur=size(all_entries(i).av_baseline,2);
            
        min_baseline_dur=min(min_baseline_dur, baseline_dur);
        min_resp_dur=min(min_resp_dur, resp_dur);  
    
        if i==1                  
            dir_vector=fly_av_vals.dir_vector;
            ndir=length(fly_av_vals.dir_vector);
        end
        i=i+1;
        clear fly_av_vals;       
    end
   
    nfly_files=i-1;
    if nfly_files==0
        return;
    end
    %compute average responses per animal
    resp_vals=zeros(nfly_files,ndir,min_resp_dur);
    bl_vals=zeros(nfly_files,ndir,min_baseline_dur);
    for j=1:nfly_files
        resp_vals(j,:,:)=all_entries(j).av_resp(:,1:min_resp_dur);
        bl_vals(j,:,:)=all_entries(j).av_baseline(:,1:min_baseline_dur);
    end
    
    average_vals.raw_resp=resp_vals;
    average_vals.raw_bl=bl_vals;
    average_vals.av_resp=squeeze(mean(resp_vals,1));
    average_vals.std_resp=squeeze(std(resp_vals,[],1));
    average_vals.av_bl=squeeze(mean(bl_vals,1));
    average_vals.std_bl=squeeze(std(bl_vals,[],1));
    average_vals.n_flies=nfly_files;      
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
