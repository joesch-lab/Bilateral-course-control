function save_all_reps_scanning_rect_raw_data_folder(curfolder)
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
        
    %% iterate thru stim files, get raw RF values
    nfiles=length(stim_files);
    ncount=0;
    logiles={};
    prfiles={};
    repi=1;
    for i=1:nfiles

        fulllog_file=fullfile(curfolder,stim_files{i});
        stimparam_i= read_scanning_rect_with_pause(fulllog_file);
        if stimparam_i.rheight~=30, continue; end

        rawdata_i=scanning_rect_one_recording_seprep(fulllog_file);
        if isempty(rawdata_i)
            continue;
        end

        %save data per rep
        for ri=1:rawdata_i.stim_info.nrep

            clear rep_data;
            rep_data.folder=curfolder;
            rep_data.logfiles=stim_files{i};
            rep_data.prfiles=rawdata_i.prname;
            rep_data.nrep_i=repi;
            rep_data.stim_info.scan_h_x=rawdata_i.stim_info.scan_h_x;
            rep_data.stim_info.scan_h_y=rawdata_i.stim_info.scan_h_y;
            rep_data.stim_info.scan_v_x=rawdata_i.stim_info.scan_v_x;
            rep_data.stim_info.scan_v_y=rawdata_i.stim_info.scan_v_y;

            rep_data.trace.htrace=rawdata_i.RF(ri).hor;
            rep_data.trace.vtrace=rawdata_i.RF(ri).ver;
            rep_data.trace.h_left_t=rawdata_i.RF(ri).hor_left_start;
            rep_data.trace.h_right_t=rawdata_i.RF(ri).hor_right_start;
            rep_data.trace.v_down_t=rawdata_i.RF(ri).ver_down_start;
            rep_data.trace.v_up_t=rawdata_i.RF(ri).ver_up_start;
            rep_data.trace.v_redsignal_raw = rawdata_i.RF(ri).ver_rf;
            rep_data.trace.h_redsignal_raw = rawdata_i.RF(ri).hor_rf;
            
            rep_data.RF.RFx=rawdata_i.RF(ri).RFx;   
            rep_data.RF.RFy=rawdata_i.RF(ri).RFy;
            rep_data.RF.RFx_pos=rawdata_i.RF(ri).RFx_pos;
            rep_data.RF.RFx_neg=rawdata_i.RF(ri).RFx_neg;
            rep_data.RF.RFy_pos=rawdata_i.RF(ri).RFy_pos;
            rep_data.RF.RFy_neg=rawdata_i.RF(ri).RFy_neg;

            rep_data.strain = strain;
            rep_data.cellgroup = cellgroup;
            rep_data.celltype = celltype;
            rep_data.datestr = datestr;
            rep_data.cellid = cellid;

            %save the rep
            if isempty(cellid)
                filename_i=['scaning_rect_',strain,'_',celltype,'_',datestr,'_rep_',num2str(repi)];
            else
                filename_i=['scaning_rect_',strain,'_',celltype,'_',datestr,'_',cellid,'_rep_',num2str(repi)];
            end
            fullfile_i=fullfile(resfolder,[filename_i,'.mat']);
            save(fullfile_i,'rep_data','-v7.3');
            make_figure_RF_one_rep(rep_data,resfolder);
            repi=repi+1;
        end
    end
end