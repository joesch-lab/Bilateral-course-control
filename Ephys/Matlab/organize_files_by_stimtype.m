parent_folder='C:\DATA\Data_for_Olga\fly_ephys\FlpND\201202\copy';

%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list


%go thru the list, analize each log file with 'scanning_dot' in the title
for li =1:length(filelist)
    fname=filelist(li).name;
    if  contains(fname,'.log','IgnoreCase', true)  
        logfile_fullname=fullfile(filelist(li).folder, fname);
        switch 1    
            case contains(fname,'curtain','IgnoreCase', true) 
                stim_folder=fullfile(filelist(li).folder,'curtain_stimulus');                
            case contains(fname,'brownian','IgnoreCase', true) 
                stim_folder=fullfile(filelist(li).folder,'brownian_dots_stimulus'); 
            case contains(fname,'scanning_rect','IgnoreCase', true) 
                stim_folder=fullfile(filelist(li).folder,'scanning_rect_stimulus'); 
            case contains(fname,'S_Vel_','IgnoreCase', true) 
                stim_folder=fullfile(filelist(li).folder,'vel_tuning_stimulus'); 
            case contains(fname,'grating','IgnoreCase', true) 
                stim_folder=fullfile(filelist(li).folder,'grating_stimulus'); 
            case contains(fname,'random','IgnoreCase', true) 
                stim_folder=fullfile(filelist(li).folder,'random_dots_stimulus');             
            otherwise
                stim_folder=fullfile(filelist(li).folder,'other'); 
                disp(fname);
        end
        
        %if the file already in the stim subfolder, skip it
        [~,current_subfolder,~]=fileparts(filelist(li).folder);
        [~,new_folder,~]=fileparts(stim_folder);
        if strcmp(current_subfolder,new_folder)
            continue;
        end
        
        %if folder does not exist, create it
        if ~exist(stim_folder)
            mkdir(stim_folder);
        end

        %find the closest by data pr file
        prfilename = find_pr_file_closest_date(logfile_fullname);
        if(prfilename==-1)
            disp(['	... in',filelist(li).folder]);
            continue;
        end
        prfilename_full=fullfile(filelist(li).folder,prfilename);

        %move log and pr files into the stim folder
        movefile(logfile_fullname, stim_folder);
        movefile(prfilename_full, stim_folder);    
    end
end


        
function prfilename = find_pr_file_closest_date(logfilename)
   slashinds=strfind(logfilename,filesep);
   folder=logfilename(1:slashinds(end));
   log_info = dir(logfilename);
   log_datevec=datevec(log_info.date);
   D = dir(folder);
   min_et=3600;
   prfilename=-1;
   for k = 3:length(D) % avoid using the first ones
       currfile = D(k).name; %file name
       if contains(currfile,'.pr','IgnoreCase', true)
           et=abs(etime(datevec(D(k).date),log_datevec));
           if et<min_et
               min_et=et;
               prfilename=currfile;
           end
       end
   end
   if min_et>10
       disp('No matching pr file has been found');
       prfilename=-1;
   end      
end
