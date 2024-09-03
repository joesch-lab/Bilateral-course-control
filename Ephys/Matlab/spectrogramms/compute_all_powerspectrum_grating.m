% compute_all_powerspectrum_from_grating

parent_folder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';

%get the list of files in all folders and subfolders
allsubfolders=['**',filesep,'*_raw_ds.mat'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
vel_pattern = 'S_Vel_';
%go thru the list, analize each MAT file
for li =1:length(filelist)
    %skip velocity tuning grating stimuli
    if contains(filelist(li).name, vel_pattern)
        continue;
    end
    disp(li);
    matresfile = fullfile(filelist(li).folder,filelist(li).name);
    grating_powerspectrum(matresfile);
end
