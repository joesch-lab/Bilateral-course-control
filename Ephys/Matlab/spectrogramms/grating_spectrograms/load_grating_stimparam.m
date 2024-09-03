function stimparam = load_grating_stimparam(logfullname)
% read the file with stimuli parameters
    f=fopen(logfullname);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'grating_stimuli'))
        disp('The log file does not contain information about grating stimuli.');
        fclose(f);
        return;
    end
    stimname=tline;

    tline=fgetl(f); %start time of the recording    

    tline=fgetl(f); %end time of the recording
    tstr=tline(length('End time: ')+1:end);
    formatIn='yyyy/mm/dd HH:MM:SS:FFF';
    stimparam.time=datevec(tstr,formatIn);
    
    tline=fgetl(f); %time per grating     
    tstr=tline(length('Time per grating: ')+1:end);    
    stimparam.gr_dur=str2double(tstr); 

    tline=fgetl(f); %Flip rate: 
    
    tline=fgetl(f);
    tstr=tline(length('Repetitions: ')+1:end);    
    stimparam.nrep=str2double(tstr); 
    
    tline=fgetl(f);
    tstr=tline(length('Speed: ')+1:end);    
    stimparam.speed=str2double(tstr); 
        
    tline=fgetl(f);
    tstr=tline(length('Baseline recording: ')+1:end);     
    stimparam.baseline_dur=uint32(str2num(tstr)); 
    
    tline=fgetl(f);
    tstr=tline(length('Spatial frequency: ')+1:end);     
    stimparam.sp_freq=str2double(tstr); 

    stimparam.VSHS='';
    if contains(stimname,'HS')
        stimparam.VSHS='HS';
    end
    if contains(stimname,'VS')
        stimparam.VSHS='VS';
    end    
    fclose(f);
end