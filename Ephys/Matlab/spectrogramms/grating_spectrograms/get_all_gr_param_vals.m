%  %check all stim durations
%  parent_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';
% 
% allsubfolders=['**',filesep,'*_grating_stimuli_*.log'];
% filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
% ndur = [];

% go thru the list and collect all folders with '_grating_stimuli_' files
files2analyze={};    
ii=1;
for li =1:length(filelist) 
    filename=fullfile(filelist(li).folder,filelist(li).name);
    if contains(filename,'\HS\') && ...
            ~(contains(filename, '_Vel') || ...
              contains(filename, '_expansion_') || ...
              contains(filename, '_contrast_levels_') ||...
              contains(filename, '_split_screen_'))
        files2analyze{ii}=filename;  
        ii=ii+1;
        stimparami=load_grating_stimparam(filename);
        ndur=[ndur;stimparami.sp_freq];
    end
end
files2analyze=files2analyze';

% files2analyze={};    
% ii=1;
% for li =1:length(filelist) 
%     filename=fullfile(filelist(li).folder,filelist(li).name);
%     if contains(filename,'\VS\') && ...
%             ~(contains(filename, '_Vel') || ...
%               contains(filename, '_expansion_') || ...
%               contains(filename, '_contrast_levels_') ||...
%               contains(filename, '_split_screen_'))
%         files2analyze{ii}=filename;  
%         ii=ii+1;
%         stimparami=load_grating_stimparam(filename);
%         ndur=[ndur;stimparami.gr_dur];
%     end
% end
% files2analyze=files2analyze';
