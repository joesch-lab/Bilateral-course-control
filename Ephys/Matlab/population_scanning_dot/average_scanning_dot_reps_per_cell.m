%go thru each folder in the parent folder and average all reps of scanning dots
% per folder. Each rep will be upscaled to the max size of all reps in the
% folder before averagin. There is no normalization or spherical
% deformatiion at this step. The max number of horizontal and vertical scan
% lines during reps will be strored together with RF.

function average_scanning_dot_reps_per_cell(parent_folder)
% parent_folder = 'D:\Repository\';

%get all files with reps
filelist=  get_all_scanning_dot_reps(parent_folder, "");
nfiles = length(filelist);

rec_folders={};

for li =1:nfiles
    [rec_folders{li},~,~]=fileparts(filelist{li});
end
rec_folders=unique(rec_folders);
nfolders=length(rec_folders);

%go thru the folders and average the reps within
for fi=1:nfolders
    curfolder=rec_folders{fi};
    group_data={};
    %get reps within
    reps_in_folder=get_all_scanning_dot_reps(curfolder, "");
    for li=1:length(reps_in_folder)
        rep_file=reps_in_folder{li};
        clear rep_data;

        %get one rep results
        load(rep_file);
        group_data.recording(li)=rep_data;
    end
    %%resize all RF to max size and average
    [maxrow,maxcol] = get_RF_max_size(group_data);
    RF_x=zeros(maxrow,maxcol);
    RF_y=zeros(maxrow,maxcol);
    nrec_real=0;
    repinds=[];
    maxHlines=0;
    maxVlines= 0;
    repfiles={};
    for ri=1:length(reps_in_folder)
        RFxi=group_data.recording(ri).RF.RFx;
        if any(isnan(RFxi(:)))
            continue;
        end
        RFyi=group_data.recording(ri).RF.RFy;
        if any(isnan(RFyi(:)))
            continue;
        end
        nrec_real=nrec_real+1;
        repfiles{nrec_real}=reps_in_folder(ri);
        repinds=[repinds,group_data.recording(ri).nrep_i];
        maxHlines=max(maxHlines,size(RFxi,1));
        maxVlines= max(maxVlines,size(RFyi,2));


        RFxi=imresize(RFxi,[maxrow,maxcol],"nearest");
        RFyi=imresize(RFyi,[maxrow,maxcol],"nearest");
        RF_x=RF_x+RFxi;
        RF_y=RF_y+RFyi;
    end
    RF_x=RF_x/nrec_real;
    RF_y=RF_y/nrec_real;

    folder_data.strain = rep_data.strain;
    folder_data.cellgroup = rep_data.cellgroup;
    folder_data.celltype = rep_data.celltype;
    folder_data.datestr = rep_data.datestr;
    folder_data.cellid = rep_data.cellid;
    folder_data.repfiles=repfiles;
    folder_data.nreps=nrec_real;
    folder_data.repinds=repinds;
    folder_data.maxHlines=maxHlines;
    folder_data.maxVlines=maxVlines;
    folder_data.RF_x=RF_x;
    folder_data.RF_y=RF_y;

    
    if isempty(folder_data.cellid)
        filename_i=['scaning_rect_',folder_data.strain,'_',folder_data.celltype,'_',folder_data.datestr,'_cellav'];
    else
        filename_i=['scaning_rect_',folder_data.strain,'_',folder_data.celltype,'_',folder_data.datestr,'_',folder_data.cellid,'_cellav'];
    end
    figfile=fullfile(curfolder,[filename_i,'.fig']);
    make_figure_RF_one_cell(folder_data,figfile);
    fullfile_i=fullfile(curfolder,[filename_i,'.mat']);
    save(fullfile_i,'folder_data','-v7.3');
end
end

function filelist=  get_all_scanning_dot_reps(parent_folder, group_labels)
    filelist ={};
        
    %get the list of files in all folders and subfolders
    allsubfolders=['**',filesep,'scaning_rect_*_rep_*.mat'];
    filelist_all = dir(fullfile(parent_folder, allsubfolders));%get list of files and folders in any subfolder
    npat=size(group_labels,1);
    %go thru the list and collect all folders with '_scanning_rect_' files
    lii=1;
    for li =1:length(filelist_all)
        fullname=fullfile(filelist_all(li).folder,filelist_all(li).name);
        for patti=1:npat
            if str_contains_strarray(fullname,group_labels(patti,:))
                filelist{lii}=fullname;
                lii=lii+1;
                break;
            end
        end
    end    
end

function TF =  str_contains_strarray(strname,strarray)
    label_len=size(strarray,2);
    TF=1;
    for i=1:label_len
        TF=TF && contains(strname,strarray(i));
    end
end

function [maxrow,maxcol] = get_RF_max_size(group_data)
    nrec=numel(group_data.recording);
    maxrow=0;
    maxcol=0;
    for ri=1:nrec
        maxrow=max(maxrow,size(group_data.recording(ri).RF.RFy,1));
        maxcol=max(maxcol,size(group_data.recording(ri).RF.RFx,2));
    end
end