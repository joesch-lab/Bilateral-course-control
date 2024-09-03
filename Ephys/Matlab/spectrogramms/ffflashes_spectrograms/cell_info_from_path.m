function [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(folderpath)
%search for the entry \HS\ or \VS\ in the path
    strain = '';
    cellgroup ='';
    celltype = '';
    datestr = '';
    cellid='';

    if ~endsWith(folderpath,filesep)
        folderpath=[folderpath,filesep]
    end

    grp_len=2;
    mypatt='[\\\/][HV][2S][\\\/]';
    k = regexp(folderpath,mypatt);

    if isempty(k) || length(k)>1
        return;
    end
    cellgroup =folderpath(k+1:k+grp_len);   

    %folder level above cell group is strain name
    [~,strain,~]=fileparts(folderpath(1:k-1));
%     k=k+grp_len+1;


    
    %folder level under cell group is cell type
    slashpatt='[\\\/]';
    start_ind = k+4;
    slashinds=start_ind+regexp(folderpath(start_ind:end),slashpatt)-1;
   if isempty(regexp(folderpath(start_ind:slashinds(1)-1),'\d{6}'))
        celltype=folderpath(start_ind:slashinds(1)-1);
        datestr='';
    else
        celltype='';
        datestr=folderpath(start_ind:slashinds(1)-1);
    end
    
    %date string is the folder level below cell type
    if length(slashinds)<2 %if the last slash is absent
        return;
    end
    if isempty(datestr)
        datestr=folderpath(slashinds(1)+1:slashinds(2)-1);
    else
        cellid=folderpath(slashinds(1)+1:slashinds(2)-1);
    end

    if ~isempty(cellid)
        return;
    end

    if length(folderpath)<=slashinds(2)
        return; %no cell id level
    end
    
    cellid=folderpath(slashinds(2)+1:slashinds(3)-1);
end
