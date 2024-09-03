function process_all_grating_powerspectrum_folders(parent_folder)
    %each folder contains recording for the same cell, hence multiple
    %recordings should be averaged
    %pattern to include
    patt_incl = '\HS\';
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*_grating_stimuli_*.log'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with '_grating_stimuli_' files
    folders2analyze={}; 
    i=1;
    for li =1:length(filelist) 
        filename=fullfile(filelist(li).folder,filelist(li).name);
        if contains(filename,patt_incl) && ...
            ~(contains(filename, '_Vel') || ...
              contains(filename, '_expansion_') || ...
              contains(filename, '_contrast_levels_') ||...
              contains(filename, '_split_screen_'))
            
            folders2analyze{i}=filelist(li).folder;    
            i=i+1;
        end
    end
    folders2analyze=unique(folders2analyze');
    disp(['Total folders with the stimulus: ',num2str(length(folders2analyze))]);
    
    %analyze and average recordings within one folder  
    for li =1:length(folders2analyze)       
        disp(li);
        disp(folders2analyze(li));
        save_grating_stimuli_powerspectrum_folder(folders2analyze{li});        
    end
% end