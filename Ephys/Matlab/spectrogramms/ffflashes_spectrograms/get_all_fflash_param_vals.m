%  %check all stim durations
 parent_folder = '\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND_DB331';

allsubfolders=['**',filesep,'*_full_field_flashes_*.log'];
filelist = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
flashdur = [];
bd_ba=[];
nrep=[];
files2analyze={};    
pattern = '';%'\HS\'
% go thru the list and collect all folders with '_full_field_flashes_' files
ii=1;
for li =1:length(filelist) 
    filename=fullfile(filelist(li).folder,filelist(li).name);
    if contains(filename,pattern) 
        files2analyze{ii}=filename;  
        ii=ii+1;
        stimparami=load_ffflashes_stimparam(filename);
        flashdur=[flashdur;stimparami.flash_dur];
        bd_ba=[bd_ba;[stimparami.tbefore, stimparami.tafter]];
        nrep=[nrep;stimparami.nrep];        
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
