%collect all reports of the missing frames in one file
parent_folder = 'D:\Repository\FlpD';
missing_patt = '*_missing_frames_info.txt';
report = fullfile(parent_folder, 'all_missing_frame_info.txt');

filename_gen = fullfile(parent_folder,'**',missing_patt);
filelist=dir(filename_gen);

fr=fopen(report,'w');
for i=1:length(filelist)
    filenamei=fullfile(filelist(i).folder,filelist(i).name);
    fprintf(fr, "%s\n", filenamei);

    f=fopen(filenamei);
    tline=fgetl(f);
    while  ~feof(f) && isempty(tline)
        tline=fgetl(f);
    end
    while  ~feof(f) && ~isempty(tline)        
        fprintf(fr, "%s\n", tline);
        tline=fgetl(f);
    end
    fclose(f);
    fprintf(fr, "%s\n", tline);
end
fclose(fr);