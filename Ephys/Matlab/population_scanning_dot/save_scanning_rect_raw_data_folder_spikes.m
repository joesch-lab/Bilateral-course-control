function save_scanning_rect_raw_data_folder_spikes(curfolder)
% curfolder='\\fs.ist.ac.at\dfsgroup\joeschgrp\Vika\EPhys\shakB_project\Repository\FlpND\VS\VS1_4\210128';
    %% get strain, cell type and date from the folder path
    [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(curfolder);
    
    %% get all stimulus files
    filelist = dir(curfolder);%get list of files
    filelist = filelist(~[filelist.isdir]);  %remove folders from list

    stim_name='_scanning_rect_';
    
    %go thru the list and collect files with'_scanning_rect_' 
    stim_files={};
    i=1;
    for li =1:length(filelist)
        match_ind = regexp(filelist(li).name, stim_name);
        if ~isempty(match_ind) &&  contains(filelist(li).name,'.log','IgnoreCase', true)  
            stim_files{i}=filelist(li).name;            
            i=i+1;
        end
    end
    
    resfolder=fullfile(curfolder,'res');
    if ~exist(resfolder,'dir')
        mkdir(resfolder);
    end
    
    total_nrep=0;
    %% iterate thru stim files, get raw RF values
    nfiles=length(stim_files);
    ncount=0;
    logiles={};
    prfiles={};
    for i=1:nfiles       
       
        fulllog_file=fullfile(curfolder,stim_files{i});
        stimparam_i= read_scanning_rect_with_pause(fulllog_file);
        if stimparam_i.rheight~=30, continue; end
                
        rawdata_i=scanning_rect_one_recording_spikes(fulllog_file); 
        if isempty(rawdata_i)            
            continue;
        end
        ncount=ncount+1;
        %structure to pool all data        
        logiles{ncount}=stim_files{i};        
        prfiles{ncount}=rawdata_i.prname;
        nreps=stimparam_i.nrep;
        total_nrep=total_nrep+nreps;
        RF(ncount).RFx = rawdata_i.RF.RFx;
        RF(ncount).RFy = rawdata_i.RF.RFy;
        RF(ncount).RFx_pos = rawdata_i.RF.RFx_pos;
        RF(ncount).RFx_neg = rawdata_i.RF.RFx_neg;
        RF(ncount).RFy_pos = rawdata_i.RF.RFy_pos;
        RF(ncount).RFy_neg = rawdata_i.RF.RFy_neg;        
    end
    
    %% average RF per folder
    % upscale RF before averaging   
    [nrow,ncol]=size(RF(1).RFx);
    RFx=zeros(nrow,ncol);
    RFy=zeros(nrow,ncol);
    RFx_pos=zeros(nrow,ncol);
    RFx_neg=zeros(nrow,ncol);
    RFy_pos=zeros(nrow,ncol);
    RFy_neg=zeros(nrow,ncol);
    for i=1:ncount
        [RFx_hres, RFy_hres]= upscale_scanning_dimensions_fullscreen(RF(i).RFx,RF(i).RFy);
        RFx=RFx+RFx_hres;
        RFy=RFy+RFy_hres;
        [RFx_hres, RFy_hres]= upscale_scanning_dimensions_fullscreen(RF(i).RFx_pos,RF(i).RFy_pos);
        RFx_pos=RFx_pos+RFx_hres;
        RFy_pos=RFy_pos+RFy_hres;
        [RFx_hres, RFy_hres]= upscale_scanning_dimensions_fullscreen(RF(i).RFx_neg,RF(i).RFy_neg);        
        RFx_neg=RFx_neg+RFx_hres;        
        RFy_neg=RFy_neg+RFy_hres;
    end
    RFx=RFx./ncount;    RFy=RFy./ncount;
    RFx_pos=RFx_pos./ncount;    RFx_neg=RFx_neg./ncount;
    RFy_pos=RFy_pos./ncount;    RFy_neg=RFy_neg./ncount;
      
    %create data structure
    clear exp_data;
    exp_data.folder=curfolder;
    exp_data.logfiles=logiles;
    exp_data.prfiles=prfiles;
    
    
    exp_data.nrep= total_nrep;
    exp_data.RF.all_RF = RF;   
    exp_data.RF.av_RFx=RFx;   exp_data.RF.av_RFy=RFy;
    exp_data.RF.av_RFx_pos=RFx_pos;   exp_data.RF.av_RFx_neg=RFx_neg;
    exp_data.RF.av_RFy_pos=RFy_pos;   exp_data.RF.av_RFy_neg=RFy_neg;  
    
    exp_data.strain = strain;
    exp_data.cellgroup = cellgroup;
    exp_data.celltype = celltype;
    exp_data.datestr = datestr;
    exp_data.cellid = cellid;       

    %save    
    if isempty(cellid)
        filename_i=['scaning_rect_spikes_',strain,'_',celltype,'_',datestr];        
    else
        filename_i=['scaning_rect_spikes_',strain,'_',celltype,'_',datestr,'_',cellid];        
    end    
    fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
    save(fullfile_i,'exp_data','-v7.3');  
end