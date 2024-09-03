function all_clevels = all_contast_levels(parent_folder)
    %collect all possible contrast levels info in the folder
    all_clevels=[];

    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*_stimuli_contrast_levels_*log.mat'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list, read params
    for li =1:length(filelist)    
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        stimparam= load(logfile_fullname);        
        all_clevels=[all_clevels; stimparam.contrast_levels];    
    end
    

    
    %go thru the list, read params
    for li =1:length(filelist)    
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        stimparam= load(logfile_fullname);        
        disp(logfile_fullname);
        disp(unique(stimparam.allconfigs(:,2)));    
    end
    
 %go thru the list, read params
    for li =1:length(filelist)    
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        stimparam= load(logfile_fullname);        
        disp(logfile_fullname);
        disp(stimparam.withgray_yn);    
    end
   

    all_clevels=unique(speed_array,'rows');
    %sort by contrast
    all_clevels_diff=all_clevels;
    all_clevels_diff(:,3)=abs(all_clevels_diff(:,1)-all_clevels_diff(:,2));
    all_clevels=sortrows(all_clevels,3,"descend");
    all_clevels=all_clevels(:,1:2);
end