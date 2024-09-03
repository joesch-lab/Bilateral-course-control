function process_all_scanning_rect_powerspectrum_folders(parent_folder)
    %each folder contains recording for the same cell, hence multiple
    %recordings should be averaged
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*_scanning_rect_*.log'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with '_scanning_rect_' files
    folders2analyze={};    
    for li =1:length(filelist)          
        folders2analyze{li}=filelist(li).folder;                           
    end
    folders2analyze=unique(folders2analyze);
    disp(['Total folders with the stimulus: ',num2str(length(folders2analyze))]);
    
    %analyze the folder, combining recordings with the same parameters    
    for li =1:length(folders2analyze)       
        disp(li);
        disp(folders2analyze(li));
        save_scanning_rect_powerspectrum_all_sides_folder(folders2analyze{li});        
    end
end