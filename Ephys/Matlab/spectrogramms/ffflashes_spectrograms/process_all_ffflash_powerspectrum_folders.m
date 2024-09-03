function process_all_ffflash_powerspectrum_folders(parent_folder)
    %each folder contains recording for the same cell, hence multiple
    %recordings should be averaged
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*_full_field_flashes_*.log'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with 'full_field_flashes' files
    folders2analyze={}; 
    i=1;
    for li =1:length(filelist)            
        folders2analyze{i}=filelist(li).folder;    
        i=i+1;
    end
    folders2analyze=unique(folders2analyze');
    disp(['Total folders with the stimulus: ',num2str(length(folders2analyze))]);
    
    %analyze and average recordings within one folder  
    for li =1:length(folders2analyze)       
        disp(li);
        disp(folders2analyze(li));
        save_ffflashes_stimuli_powerspectrum_folder(folders2analyze{li});
        save_ffflashes_mean_resps_folder(folders2analyze{li});
    end
end