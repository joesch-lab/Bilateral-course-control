function [strain,cellgroup,celltype,datestr,cellid] = cell_info_from_path(folderpath)
%search for the entry \HS\ or \VS\ in the path
    strain = '';
    cellgroup ='';
    celltype = '';
    datestr = '';
    cellid='';
    grp_len=2;
    mypatt='[\\\/][HV][2S][\\\/]';
    k = regexp(folderpath,mypatt);

    if isempty(k) || length(k)>1
        return;
    end
    cellgroup =folderpath(k+1:k+grp_len);
    
    %foder level above cell group is strain name
    [~,strain,~]=fileparts(folderpath(1:k-1));
    
    %foder level under cell group is cell type
    slashpatt='[\\\/]';
    start_ind = k+4;
    slashinds=start_ind+regexp(folderpath(start_ind:end),slashpatt)-1;
    celltype=folderpath(start_ind:slashinds(1)-1);
    
    %date string is the folder level below cell type
    if length(slashinds)<2 %if the last slash is absent
        slashinds(2)=length(folderpath)+1;
    end
    datestr=folderpath(slashinds(1)+1:slashinds(2)-1);
    
    if length(folderpath)<slashinds(2)+1
        return; %no cell id level
    end
    
    if length(slashinds)<3 %if the last slash is absent
        slashinds(3)=length(folderpath)+1;
    end
    cellid=folderpath(slashinds(2)+1:slashinds(3)-1);
end
