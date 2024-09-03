curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\';

subflds = fullfile(curfolder,'\**\vel_tuning_*.mat');

filelist = dir(subflds);%get list of files
   

    for li =1:length(filelist)
        fulllogname=fullfile(filelist(li).folder,filelist(li).name);
        disp(fulllogname);
        delete(fulllogname);
    end