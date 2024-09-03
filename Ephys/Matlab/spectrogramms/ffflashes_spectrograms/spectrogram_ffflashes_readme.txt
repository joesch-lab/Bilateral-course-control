Short description of the files in the folder:

process_all_ffflash_powerspectrum_folders.m
this is the main script. It collects all the folders that contain _full_filed_flashes_ recordings. 
It launches the analysis per folder. The results are averaged per folder.

save_ffflashes_stimuli_powerspectrum_folder.m
for a given folder (aka animal), compute power spectrum ananalysis per rep and average.

trace_on_of_flashes_one_recording.m 
for a given log file, finds the corresponding *.pr recording and returns 
the raw traces during the on-flash and period after the flash.               

make_group_plots.m    
collecte all '*flashes_on_off_pwspectrum_*.mat' files with the results of powerspectrum analysis, 
group them within the specific groups (e.g strain types) and make plots

load_ffflashes_stimparam.m      
reads individual flash log file and returns params

cell_info_from_path.m     
given a file path, extracts strain, cell group, cell type, date, and cell_id values

get_all_ffflashes_param_vals.m
for all fullfield flash files saves the values of a specified param along with the name of the log file.
Use this script to check all possible values of a specific parameter.  
  