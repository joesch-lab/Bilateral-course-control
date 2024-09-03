function stimparam = read_velocity_tuning_stimuli_log(logfile)
    stimparam={};
    % read the file with stimuli parameters
    f=fopen(logfile);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'_Vel'))
        disp('The log file does not contain information about velocity selectivity stimuli.');
        fclose(f);
        return;
    end
    
    if isempty(strfind(tline,'HS'))
        isHS=0;
    else
        isHS=1;
    end
    
    tline=fgetl(f);
    tstr=tline(length('Start time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    tstart=datevec(tstr,formatIn);

    tline=fgetl(f);
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    tend=datevec(tstr,formatIn);
    
    tline=fgetl(f);
    t_str=tline(length('Time per grating: ')+1:end);     
    tmove=uint32(str2double(t_str)); 
    
    tline=fgetl(f);
    ifistr=tline(length('Flip rate: ')+1:end);    
    ifi=str2double(ifistr);
    
    tline=fgetl(f);
    nrepstr=tline(length('Repetitions: ')+1:end);    
    nrep=str2num(nrepstr); 

    speed_array=[];    
    tline=fgetl(f);
    while ~isempty(strfind(tline,'Speed'))
        tstr=tline(length('Speed: ')+1:end);    
        speedi=str2double(tstr); 
        speed_array=[speed_array,speedi];
         tline=fgetl(f);
    end
        
    %Baseline recording: skipped
    
    tline=fgetl(f);
    sf_str=tline(length('Spatial frequency: ')+1:end);     
    sp_freq=str2double(sf_str); 
    
    speed_indices=[];
    tline=fgetl(f);
    while ~isempty(strfind(tline,'Random_Number'))
        tstr=tline(length('Random_Number: ')+1:end);    
        speedi=str2num(tstr); 
        speed_indices=[speed_indices,speedi];
         tline=fgetl(f);
    end
    
    fclose(f);
    stimparam.nrep=nrep;
    stimparam.tmove=tmove;
    stimparam.speed_array=speed_array;
    stimparam.speed_indices=speed_indices;
    stimparam.sp_freq=sp_freq;
    stimparam.ifi=ifi;
    stimparam.isHS=isHS;    
end