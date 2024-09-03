function stimparam = load_ffflashes_stimparam(logfullname)
% read the file with stimuli parameters
    f=fopen(logfullname);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end

    if isempty(strfind(tline,'full_field_flashes'))
        disp('The log file does not contain information about full_field_flashes stimuli.');
        fclose(f);
        return;
    end
    stimname=tline(length('Stimuli: ')+1:end);

    tline=fgetl(f); %start time of the recording  
    tline=fgetl(f); %end time of the recording
    tline=fgetl(f); %Flip rate: 

    tline=fgetl(f);
    tstr=tline(length('Number of reps: ')+1:end); 
    stimparam.nrep=str2double(tstr); 
    
    tline=fgetl(f); %time per flash     
    tstr=tline(length('tflash: ')+1:end);    
    stimparam.flash_dur=str2double(tstr); 

    tline=fgetl(f); %background before flash     
    tstr=tline(length('tbefore: ')+1:end);    
    stimparam.tbefore=str2double(tstr); 

    tline=fgetl(f); %background after flash   
    tstr=tline(length('tafter: ')+1:end);    
    stimparam.tafter=str2double(tstr); 

    tline=fgetl(f); %flash color  
    tstr=tline(length('Flash_color: [')+1:end-1);
    tstr=strsplit(tstr);    
    stimparam.flash_color=cellfun(@str2num,tstr); 

    tline=fgetl(f); %background color  
    tstr=tline(length('BG_color: [')+1:end-1);
    tstr=strsplit(tstr);    
    stimparam.bg_color=cellfun(@str2num,tstr); 

    fclose(f);
end