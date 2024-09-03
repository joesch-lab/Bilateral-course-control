function [speed_array, sp_freq_array] = all_speed_spfreq_param_vel_tuning(parent_folder)
    %collect all possible speed and sp frequency info
    speed_array=[];
    sp_freq_array=[];

    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'*S_Vel_*.log'];
    filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    
    %go thru the list, read params
    for li =1:length(filelist)    
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        stimparam= read_velocity_tuning_stimuli_log(logfile_fullname);
        sp_freq_array=[sp_freq_array,stimparam.sp_freq];
        speed_array=[speed_array,stimparam.speed_array];        
    end
    
    speed_array=unique(speed_array);
    sp_freq_array=unique(sp_freq_array);
end