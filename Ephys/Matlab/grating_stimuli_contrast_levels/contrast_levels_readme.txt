Short description of the files in the folder:

process_all_contrast_levels_raw_data.m              
this is the main script. It collects all the folders that contain grating of various contrasts recordings. 
It launches the analysis of per folder. The results are stored per folder and per different 
set of parameters.

save_contrast_levels_raw_data_folder.m
for a given folder (aka animal), analyze and group the grating with different contrasts
recordings with the same parameters. The results of all stimuli files 
from the folder are collected and averaged. The results are stored in the 
mat file per folder and per set of parameters.

contrast_levels_one_recording.m 
analyze one recording of grating with different contrast (1 .log+ 1 .pr)
results (full traces and baselines) are returnes. Nothing is saved on the disk.               

analyze_grating_stimuli_contrast_levels_res_per_file.m  
This is the original script for running thru all grating+contrasts files in 
the folder and subfolders.
This script will analyze each recording individually and produce plots for each log files.

all_contast_levels.m                                
read all grating+contrasts log files in the folder (e.g. whole repository)
and analyse sets of parameters

make_group_plots_contrast_levels.m  
collecte all grating_contrast_levels*.mat files with the results of grating+contrasts analysis, 
group them within the specific groups (e.g strain types) and make plots

cell_info_from_path.m     
given a file path, extracts strain, cell group, cell type, date, and cell_id values
     
clean_old_res.m                         
collects all mat files starting with grating_contrast_levels*.mat and deletes them

grating_stimuli_contrast_levels.m                   
stimulus file

openpr.m                  
reads pr file with raw voltage data   