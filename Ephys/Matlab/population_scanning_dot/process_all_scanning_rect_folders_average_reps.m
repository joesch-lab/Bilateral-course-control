function process_all_scanning_rect_folders_average_reps(parent_folder)
    %each folder contains recording for the same cell, hence multiple
    %recordings should be averaged
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'scaning_rect_*_rep_*.mat'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with '_scanning_rect_' files
    folders2analyze={};    
    for li =1:length(filelist)          
        folders2analyze{li}=filelist(li).folder;                           
    end
    folders2analyze=unique(folders2analyze);
    disp(['Total folders with the reps: ',num2str(length(folders2analyze))]);
    
    %analyze the folder, averaging reps for the same cell   
    for li =1:length(folders2analyze)       
        disp(li);
        disp(folders2analyze(li));
        average_scanning_dot_reps_per_cell_w_directions(folders2analyze{li}, 0);
        average_scanning_dot_reps_per_cell_w_directions(folders2analyze{li}, 1);
    end
end