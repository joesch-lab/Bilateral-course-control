Short description of the files in the folder:

process_all_scanning_rect_powerspectrum_folders.m
this is the main script. It collects all the folders that contain scanning rectangle recordings. 
It launches the analysis per folder. The results are averaged per folder.

save_scanning_rect_powerspectrum_folder.m
for a given folder (aka animal), compute power spectrum ananalysis per rep and average.

trace_scanning_rect_one_recording_left_part.m 
returns the concatenated trace for those parts of the scanning rect stimulus, when the rect
was in the left third of the screen(1 .log+ 1 .pr)               

make_group_plots.m    
collecte all scaning_rect_pwspectrum_*.mat files with the results of powerspectrum analysis, 
group them within the specific groups (e.g strain types) and make plots

read_scanning_rect_with_pause.m      
reads individual scanning rect log file and returns params

cell_info_from_path.m     
given a file path, extracts strain, cell group, cell type, date, and cell_id values

save_anystim_powerspectrum_folder.m
for the log file specified with 'stim_root' (eg stimuli_full_field_flashes) 
in the current folder, find the matching pr file and compute power spectrum analysis
for the entire recording independetly on the parameters of the stimulus.
This script will be used primary for analysing shakB2 recordings
  
   