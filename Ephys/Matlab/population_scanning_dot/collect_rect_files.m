    parent_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';
    %each folder contains recording for the same cell, hence multiple
    %recordings should be averaged
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*scanning_rect*rep*.mat'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list and collect all folders with '_scanning_rect_' files
    folders2analyze={};    
    for li =1:length(filelist)          
        folders2analyze{li}=filelist(li).folder;                           
    end
    folders2analyze=unique(folders2analyze);
    disp(['Total folders with the stimulus: ',num2str(length(folders2analyze))]);
    
    