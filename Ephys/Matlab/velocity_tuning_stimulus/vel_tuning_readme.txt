Short description of the files in the folder:

process_all_vel_tuning_raw_data.m
this is the main script. It collects all the folders that contain velocity tuning recordings. 
It launches the analysis of per folder. The results are stored per folder and per different 
set of parameters.

save_velocity_tuning_raw_data_folder.m
for a given folder (aka animal), analyze and group the vel tuning recordings with the same parameters.

vel_tuning_one_recording.m 
analyze one velocity tuning recording (1 .log+ 1 .pr)               

analyze_vel_tuning_stimulus.m
This is the original script for running thru all velocity tuning log files in the folder and subfolders.
This script will analyze each recording individually and produce plots for each log files.

all_speed_spfreq_param_vel_tuning.m
read all velocity tuning log files in the folder (e.g. whole repository)
returns the array of all speeds and spatial frequences

make_group_plots.m    
collecte all vel_tuning_*.mat files with the results of velocity tuning analysis, 
group them within the specific groups (e.g strain types) and make plots

read_velocity_tuning_stimuli_log.m      
reads individual velocity tuning log file and returns params

cell_info_from_path.m     
given a file path, extracts strain, cell group, cell type, date, and cell_id values
     
clean_old_res.m                         
collects all mat files starting with vel_tuning_* and deletes them

   