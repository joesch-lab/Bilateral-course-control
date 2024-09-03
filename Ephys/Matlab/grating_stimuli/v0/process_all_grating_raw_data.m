function process_all_grating_raw_data(parent_folder)

    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*.*'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  %remove folders from list


    %go thru the list, analize each log file with 'grating_stimuli' in the title
    stim_name='grating_stimuli_[\d{4}]';
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'.log','IgnoreCase', true)  
            logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
            disp(logfile_fullname);
            save_grating_raw_data(logfile_fullname);
        end
    end
end