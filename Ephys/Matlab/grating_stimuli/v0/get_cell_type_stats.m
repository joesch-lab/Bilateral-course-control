%collect all the processing results for the given cell type, speed and
%spatial frequency

%for a given cell type print how many files of each strain, for each speed and spatial frequency
function get_cell_type_stats(parent_folder,celltype)
    %collect all the files with the results 
     %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*.*'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    %go thru the list, analize each log file with 'scanning_dot' in the title
    stim_name='grating_stimuli_[\d{4}]';    
    all_entries=array2table(zeros(0,4), 'VariableNames',{'strain_type','speed','sp_freq','count'});
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'_raw_ds.mat','IgnoreCase', true)  
            matfile=[filelist(li).folder, filesep,filelist(li).name];
            %load each file and check the file type
            load(matfile);
            %now the structure exp_data is available
            if ~contains(exp_data.cell_type, celltype, 'IgnoreCase', true) 
                continue;
            end
            %for matching files collect speed and spatial frequency stats
            all_entries= [all_entries; {exp_data.strain_type, exp_data.speed, exp_data.sp_freq, 1}];
            
            clear exp_data;
        end
    end
    
    exp_summary = groupsummary(all_entries,{'strain_type','speed','sp_freq'})
    mat_file = [celltype,'_all_conditions.mat'];
    mat_file =fullfile(parent_folder,mat_file);
    save(mat_file,'exp_summary');    
end