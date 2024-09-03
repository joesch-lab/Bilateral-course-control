Short description of the files in the folder:

process_all_grating_powerspectrum_folders.m
this is the main script. It collects all the folders that contain _grating_stimuli_ recordings. 
It launches the analysis per folder. The results are averaged per folder.

save_grating_stimuli_powerspectrum_folder.m
for a given folder (aka animal), compute power spectrum ananalysis per rep and average.

trace_pd_nd_grating_one_recording.m 
for a given log file, finds the corresponding *.pr recording and returns 
the raw traces during the pd/nd motion and the baseline.               

make_group_plots.m    
collecte all '*grating_pd_nd_pwspectrum_*.mat' files with the results of powerspectrum analysis, 
group them within the specific groups (e.g strain types) and make plots

load_grating_stimparam.m      
reads individual grating log file and returns params

cell_info_from_path.m     
given a file path, extracts strain, cell group, cell type, date, and cell_id values

get_all_gr_param_vals.m
for all grating files saves the values of a specified param along with the name of the log file.
Use this script to check all possible values of a specific parameter.  
  