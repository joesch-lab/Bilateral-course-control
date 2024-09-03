% function make_population_plot(parent_folder)
%iterates thru subfolders and collects the values from the mat files. If
%there are more than one matfile in the folder, values from several mat
%files will be averaged

parent_folder='C:\DATA\Data_for_Olga\fly_ephys\pop_plot';

allsubfolders=['**',filesep,'*.*'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

traces={};


%go thru the list, analize each log file with 'scanning_dot' in the title
for li =1:length(filelist)
    if  contains(filelist(li).name,'scanning','IgnoreCase', true) && contains(filelist(li).name,'.log','IgnoreCase', true)  
        logfile_fullname=[filelist(li).folder, filesep,filelist(li).name];
        %find the closest pr file
        prfilename = find_pr_file_closest_date(logfile_fullname);
        folder=filelist(li).folder;

    end
end
% end