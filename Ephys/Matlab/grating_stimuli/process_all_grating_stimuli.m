function process_all_grating_stimuli(parent_folder)
    %each folder contains recording for the same cell, hence multiple
    %recordings should be averaged
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*stimuli_grating_stimuli_*.log'];
    stim_name='grating_stimuli_[\d{4}]';  

    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with'grating_stimuli_' files
    folders2analyze={};
    i=1;
    for li =1:length(filelist)    
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) 
            folders2analyze{i}=filelist(li).folder;            
            i=i+1;  
        end
    end
    folders2analyze=unique(folders2analyze);
    disp(['Total folders with the stimulus: ',num2str(length(folders2analyze))]);
    
    %analyze the folder, combining recordings with the same parameters    
    for li =1:length(folders2analyze)  
        disp(li);
        disp(folders2analyze(li));
        save_grating_stimuli_raw_data_folder(folders2analyze{li});        
    end
end