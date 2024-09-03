function params = read_scanning_rect_with_pause(logfile)

    % read the file with stimuli parameters
    f=fopen(logfile);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'scanning_rect'))
        disp('The file does not contain information about scanning rect stimulus');
        fclose(f);
        return;
    end

%     tline=fgetl(f); %Start time:
%     tline=fgetl(f); %End time:
    params = struct('tstart',-1,'tend',-1,'ifi',-1,'wall_offset',-1,'nrep',-1,'rwidth',-1,'rheight',-1,'dx',-1,'freezefr',-1,'tbefore',-1,'tafter',-1,'nfr',-1);
         
    while ~feof(f)
        tline=fgetl(f);
        params = init_param(tline, params);
    end  
    
    fclose(f);
end

function params = init_param(strline, params)
    semicolon = strfind(strline,':')+1;
    varname=strline(1:semicolon);
    varval=strline(semicolon+1:end);
    
    switch varname
        case 'Start time: '
            params.tstart=datevec(varval,'YYYY/mm/dd HH:MM:SS:FFF'); %start time
        case 'End time: '
            params.tend=datevec(varval,'YYYY/mm/dd HH:MM:SS:FFF'); %start time
        case 'Flip rate: '
            params.ifi=str2double(varval); %flip rate : frames per second
        case 'Wall offset: '
            nn = strsplit(varval, {' '});
            nn= nn(~cellfun('isempty',nn));  
            params.wall_offset =cellfun(@str2double, nn)'; %wall offset
        case 'Number of repetitions: '
            params.nrep =str2double(varval); 
        case 'Rectangle width: '
            params.rwidth = str2double(varval);
        case 'Rectangle height: '
            params.rheight = str2double(varval);
        case 'Speed: '
            params.dx = str2double(varval);
        case 'Number of frames to freeze: '
            params.freezefr = uint32(str2double(varval));
        case 'Background time before: '
            params.tbefore = str2double(varval);
        case 'Background time after: '
            params.tafter = str2double(varval);
        case 'Frames presented: '
            params.nfr = str2double(varval); 
    end     
end