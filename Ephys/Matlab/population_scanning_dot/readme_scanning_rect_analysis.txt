Short description of the files in the folder:

process_all_scanning_rect_folders.m
this is the main script. It collects all the folders that contain scanning rectangle recordings. 
It launches the analysis of per folder. The results are stored per folder. If there are 
several reps or recordings per folder they will be averaged without normalization 
and the resulting RF saved in the res subfolder

save_scanning_rect_raw_data_folder.m
for a given folder (aka animal), analyze and average all scanning rectangle recordings.

scanning_rect_one_recording.m 
analyze one scanning rectangle recording (1 .log+ 1 .pr)               

files used by scanning_rect_one_recording.m:
     - openpr_flatten_translate.m
            opens the pr file and flattens the mean in case of the drift

     - reconstruct_scanning_rect_with_pause_fullres.m
        reconstructs the array of stimulus frames (x,y,dx,dy) per frame

    - read_scanning_rect_with_pause.m
        reads the log file and returns stimulus parameters

files used by save_scanning_rect_raw_data_folder.m:
    - upscale_scanning_dimensions_fullscreen.m
        before averaging resutls from different recordings, 
        upscale the results of one recording to the size of the full screen
        (one recording returns the average voltage values per scanning 
        lines of the stimulus, hence need to expand scanning lines to the full screen)

    - cell_info_from_path.m     
        given a file path, extracts strain, cell group, cell type, date, and cell_id values
     
different plotting tools:

    - scanning_rect_group_plots.m    
        collecte all scaning_rect_*.mat files with the results of scanning rect analysis, 
        group them within the specific groups (e.g strain types) and make plots
     - make_deformed_RF_figures.m             
        deforms the full size raw RF and makes plots by essentiall calling the following 3 functions:
        deform_plot_uniform_sampling, makeRFfigure, makeRFfigure_with_histograms
     - deform_plot_uniform_sampling.m
        deforms the full size raw RF 
    - makeRFfigure.m
        makes a RF figure: optic flow + heat map
    - makeRFfigure_with_histograms.m         
        makes a RF figure: optic flow + histograms of average arrow length
        across horizontal and vertical directions
    - makeRF_indiv_compfigure.m
        figure of heatmaps of separate components: voltage responses 
        during positive and negative directions of movements in X/Y directions
    - make_plots_from_folder_data.m          
        loads the data from the file with the resutls for the stimulus per animal
        (e.g. after executing save_scanning_rect_raw_data_folder.m). 
        Plots RF for this specific recording.

USE CASES:
CASE 1: If you want to check the RF of one recording. Make sure that the 
        .log and .pr files are in one folder, e.g. 'C:\animal1'.
        Your actions:
        Step 1. Run the script executing save_scanning_rect_raw_data_folder('C:\animal1');
        This will produce the file 'scaning_rect_*.mat' in the 
        folder 'C:\animal1\res'. 
        Step 2. Open the script make_plots_from_folder_data.m and adjust the 
        values of the full path of the matfile with the results ('scaning_rect_*.mat)  
        as well as the folder where the RF figures will be stored. Run the script.
CASE 2: You want to process many experiments and make group plots, 
        such as average RF for a specific strain and a specific cell type.
        Your actions:
        Step 1. define the root folder, where all the experiments, 
        which you want to include in the analysis, are (i.e. 
        parent_folder = 'C:\Repository', or parent_folder = 'C:\Repository\FlpD'). 
        Run the script process_all_scanning_rect_folders(parent_folder).
        Step 2. Open the script scanning_rect_group_plots.m and update 
        the values of the variables:
         - root folder where to search for the mat files with the results 
        per animal from the previous step, i.e. population_folder = 'C:\Repository\';
        - folder where to store figures, i.e. resfolder = 'C:\Repository\scanning_rect_res';
        - groups on how to group the data when making average plots, i.e.
        group_labels=[["HSE","FlpD"];["HSE","FlpND"]]; (this defines two groups: 
        average RF for HSE cells of the strain FlpD and average 
        RF for HSE cells for the strain FlpND)




Archived old style of the analysis:
scanning_dot_one_recording_oldstyle.m  
              
                         


   