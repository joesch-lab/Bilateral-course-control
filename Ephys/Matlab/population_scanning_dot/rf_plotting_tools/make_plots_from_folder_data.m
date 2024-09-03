%use this script to load full res data from one folder and make plots with the given number of bins
function make_plots_from_folder_data(parent_folder)
    filelist_all = dir(fullfile(parent_folder, ['**',filesep,'res', filesep,'scaning_rect_*_','*.mat']));%get list of files and folders in any subfolder
    nfiles = length(filelist_all);

    for fi=1:nfiles
        datamatfile=fullfile(filelist_all(fi).folder, filelist_all(fi).name)
        load(datamatfile);
        titlestr=[exp_data.strain,' ', exp_data.celltype, ' ', exp_data.datestr, ' ', exp_data.cellid];
        [~,figbase,~]=fileparts(datamatfile);
        
        figbase=fullfile(filelist_all(fi).folder,figbase);
        make_deformed_RF_figures(exp_data.RF.av_RFx,exp_data.RF.av_RFx, 15, titlestr,figbase);
        makeRF_indiv_compfigure(exp_data.RF.av_RFx_pos, exp_data.RF.av_RFx_neg, exp_data.RF.av_RFy_pos, exp_data.RF.av_RFy_pos, 15, titlestr, figbase);
        close all;
    end
