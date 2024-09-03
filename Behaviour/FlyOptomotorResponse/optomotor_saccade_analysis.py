import os
import numpy as np
import matplotlib.pyplot as plt
from smooth_or_saccadic_optomotor_response import *
from remove_bad_trials import *
from trajectory_straightness_analysis import *
from reject_nonresponsive_trials import *
from draw_trace_plot import *

def analyse_saccade(Dir, filename=None, sex='', p=3, condn_indices=[2,1,0], ignore_values = [], variables=[]):
    orientation = None
    if filename==None:
        filename = Dir.split(r'\'')[0].split('\\')[-1]
    print(filename+'_opto_response_'+sex+'_speed'+'.npy')
    for files in os.listdir(Dir):
        if files.endswith(filename+'_opto_response_'+sex+'_speed'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response_speed = np.load(path, allow_pickle=True)
        if files.endswith(filename+'_opto_response_'+sex+'_angular_speed'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response_angular_speed = np.load(path, allow_pickle=True)
        if files.endswith(filename+'_opto_response_'+sex+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            optomotor_response = np.load(path, allow_pickle=True)
        elif files.endswith(filename+'_opto_response_'+sex+'_saccade_peak'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            saccade_peak = np.load(path, allow_pickle=True)
        elif files.endswith(filename+'_opto_response_'+sex+'_saccade_turn'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            saccade_turn = np.load(path, allow_pickle=True)
        elif files.endswith(filename+'_opto_response_'+sex+'_pos'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_pos = np.load(path, allow_pickle=True)
        elif files.endswith(filename+'_opto_response_'+sex+'_orientation'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            orientation = np.load(path, allow_pickle=True)
    
    if type(orientation) != np.ndarray:
        orientation = np.zeros(optomotor_response.shape)
    if 'monocular' in filename:
        optomotor_response, opto_response_speed, opto_response_angular_speed, saccade_peak, saccade_turn, opto_pos, orientation = \
            remove_bad_trials([optomotor_response, opto_response_speed, opto_response_angular_speed, saccade_peak, saccade_turn, opto_pos, orientation], p=p, strict=True)
    else:
        optomotor_response, opto_response_speed, opto_response_angular_speed, saccade_peak, saccade_turn, opto_pos, orientation = \
            remove_bad_trials([optomotor_response, opto_response_speed, opto_response_angular_speed, saccade_peak, saccade_turn, opto_pos, orientation], p=p)
    Dir_saccades = os.path.join(Dir, 'Saccades')
    if not os.path.isdir(Dir_saccades):
        os.mkdir(Dir_saccades)
    Dir_processed_data = os.path.join(Dir, 'ProcessedData')
    if not os.path.isdir(Dir_processed_data):
        os.mkdir(Dir_processed_data)

    straightness_list_stim, straightness_list_prestim, conditions, fwd_walking_segments_list_stim, fwd_walking_segments_list_pre = \
    straightness_analysis_hemi(saccade_peak, opto_response_speed, opto_pos, optomotor_response, opto_response_angular_speed, orientation,
        Dir_saccades, filename, p, condn_indices=condn_indices, ignore_values=ignore_values, variables=variables)
    straightness_data = {'variables':variables, 'conditions':conditions, 'stim':straightness_list_stim, 'prestim':straightness_list_prestim,
                         'walk_stim':fwd_walking_segments_list_stim, 'walk_pre':fwd_walking_segments_list_pre}
    np.save(os.path.join(Dir_processed_data, filename+'_straightness.npy'), straightness_data, allow_pickle=True)

    basic_saccade_analysis(saccade_peak, saccade_turn, optomotor_response, opto_response_speed, opto_response_angular_speed,orientation,
        Dir_saccades, filename, sex='', p=p, condn_indices=condn_indices, ignore_values=ignore_values, variables=variables)
    return 1


def basic_saccade_analysis(saccade_peak, saccade_turn, optomotor_response, opto_response_speed, opto_response_angular_speed, orientation, Dir, filename, sex='', p=3,
                           condn_indices=[2, 1 ,0], ignore_values=[], variables=[]):
    Dir_processed_data = os.path.join(os.path.dirname(Dir), 'ProcessedData')
    if not os.path.isdir(Dir_processed_data):
        os.mkdir(Dir_processed_data)
    ##
    indices, rejected_indices = reject_nonresponsive_trials(opto_response_angular_speed, opto_response_speed, p=p, stim_period=stim_period)
    optomotor_response = optomotor_response[indices]
    opto_response_speed = opto_response_speed[indices]
    opto_response_angular_speed = opto_response_angular_speed[indices]
    saccade_peak = saccade_peak[indices]
    saccade_turn = saccade_turn[indices]
    orientation = orientation[indices]
    ##
    opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=optomotor_response[:, p + 4:])
    ## add the expected response collumn after p collumns, defaut is 1
    optomotor_response, exp_type, labels = get_expected_response_hemi(optomotor_response, p=p, translation='N')
    ##
    saccade_peak[:, 0] = optomotor_response[:, p]
    # indices, conditions = indices_by_condition([exp_type], ignore_values=[[5]]) ## take FtB, BtF and full_rotation experiments
    ## multiply by expected direction such that expected values are positive and opposite responses are negative
    opto_response_list = np.multiply(opto_response_list, optomotor_response[:, p].reshape(-1, 1) * -1)
    opto_response_angular_speed_list = np.multiply(opto_response_angular_speed[:, p+3:], optomotor_response[:, p].reshape(-1, 1) * -1)
    ##
    ## cluster by BtF, FtB and full rotation and plot separtely
    data, conditions, indices = split_data_by_conditions([saccade_peak], [exp_type], ignore_values=[[5]])
    indices, conditions = indices_by_condition([exp_type], ignore_values=[[5]]) ## take FtB, BtF and full_rotation experiments
    print(conditions, labels)
    ## saccade turn distribution
    ## multiply with expected response direction to fix saccade sign
    saccade_turn_dir_corrected = np.multiply(saccade_turn, saccade_peak[:, 0].reshape(-1, 1) * -1)
    saccade_peak_dir_corrected = np.multiply(saccade_peak, saccade_peak[:, 0].reshape(-1, 1) * -1)

    syn_saccade_times, anti_saccade_times, anti_saccade_response_list, syn_saccade_response_list, pre_anti_saccade_response_list, pre_syn_saccade_response_list, \
        syn_saccade_turns, anti_saccade_turns, syn_saccade_peaks, anti_saccade_peaks, pre_syn_saccade_turns, pre_anti_saccade_turns, pre_syn_saccade_peaks, pre_anti_saccade_peaks = time_between_saccades(data, [saccade_turn_dir_corrected[index] for index in indices], p=p, start=0.0, duration=5.0) ## syn and anti are determined by the value of saccade_peak[;, 0]

    saccade_time_data = {'variables': labels, 'conditions': conditions, 'syn': saccade_timewise_syn, 'anti': saccade_timewise_anti, 'syn_saccade_dist': syn_saccade_response_list,
                         'anti_saccade_dist': anti_saccade_response_list, 'syn_pre_saccade_dist': pre_syn_saccade_response_list, 'anti_pre_saccade_dist': pre_anti_saccade_response_list,
                         'syn_saccade_dist_early': syn_saccade_response_list_early, 'anti_saccade_dist_early': anti_saccade_response_list_early, 'syn_pre_saccade_dist_early': pre_syn_saccade_response_list_early,
                         'anti_pre_saccade_dist_early': pre_anti_saccade_response_list_early, 'syn_saccade_dist_late': syn_saccade_response_list_late, 'anti_saccade_dist_late': anti_saccade_response_list_late,
                         'syn_pre_saccade_dist_late': pre_syn_saccade_response_list_late, 'anti_pre_saccade_dist_late': pre_anti_saccade_response_list_late,
                         'syn_saccade_turns_early':syn_saccade_turns_early, 'anti_saccade_turns_early':anti_saccade_turns_early, 'syn_saccade_turns_late':syn_saccade_turns_late,
                         'anti_saccade_turns_late':anti_saccade_turns_late, 'syn_saccade_turns':syn_saccade_turns, 'anti_saccade_turns':anti_saccade_turns}
    np.save(os.path.join(Dir_processed_data, filename+'_timewise_saccade_hemi.npy'), saccade_time_data, allow_pickle=True)

    saccade_time_data_per_fly = {'variables': labels, 'conditions': conditions, 'syn_saccade_dist': [], 'anti_saccade_dist': [], 'syn_pre_saccade_dist': [],'anti_pre_saccade_dist': [],
                                 'syn_saccade_turns': [], 'anti_saccade_turns': [], 'syn_saccade_peaks': [], 'anti_saccade_peaks': [],
                                 'pre_syn_saccade_turns': [], 'pre_anti_saccade_turns': [], 'pre_syn_saccade_peaks': [], 'pre_anti_saccade_peaks': [],
                                 'syn_saccade_turns_allstimuli': [], 'anti_saccade_turns_allstimuli': [], 'syn_saccade_peaks_allstimuli': [], 'anti_saccade_peaks_allstimuli': [],
                                 'pre_syn_saccade_turns_allstimuli': [], 'pre_anti_saccade_turns_allstimuli': [], 'pre_syn_saccade_peaks_allstimuli': [], 'pre_anti_saccade_peaks_allstimuli': []}
    saccade_turn_list = [saccade_turn_dir_corrected[index] for index in indices]
    fly_id = optomotor_response[:, p + 2]
    fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]

    syn_saccade_turns_per_fly_all_stimuli = []
    anti_saccade_turns_per_fly_all_stimuli = []
    syn_saccade_peaks_per_fly_all_stimuli = []
    anti_saccade_peaks_per_fly_all_stimuli = []
    pre_syn_saccade_turns_per_fly_all_stimuli = []
    pre_anti_saccade_turns_per_fly_all_stimuli = []
    pre_syn_saccade_peaks_per_fly_all_stimuli = []
    pre_anti_saccade_peaks_per_fly_all_stimuli = []
    for i in range(len(data[0])):
        syn_saccade_response_list_per_fly = []
        anti_saccade_response_list_per_fly = []
        pre_anti_saccade_response_list_per_fly = []
        pre_syn_saccade_response_list_per_fly = []
        syn_saccade_turns_per_fly = []
        anti_saccade_turns_per_fly = []
        pre_syn_saccade_turns_per_fly = []
        pre_anti_saccade_turns_per_fly = []
        syn_saccade_peaks_per_fly = []
        anti_saccade_peaks_per_fly = []
        pre_syn_saccade_peaks_per_fly = []
        pre_anti_saccade_peaks_per_fly = []
        for ids in np.unique(fly_id):
            fly_indices = np.argwhere(fly_id_indices[i] == ids).flatten()
            if fly_indices.shape[0] < 5:
                #If the number of trials is less than 5 then just put NaN
                syn_saccade_response_list_per_fly.append(np.nan)
                anti_saccade_response_list_per_fly.append(np.nan)
                pre_anti_saccade_response_list_per_fly.append(np.nan)
                pre_syn_saccade_response_list_per_fly.append(np.nan)
                syn_saccade_turns_per_fly.append(np.nan)
                anti_saccade_turns_per_fly.append(np.nan)
                pre_syn_saccade_turns_per_fly.append(np.nan)
                pre_anti_saccade_turns_per_fly.append(np.nan)
                syn_saccade_peaks_per_fly.append(np.nan)
                anti_saccade_peaks_per_fly.append(np.nan)
                pre_syn_saccade_peaks_per_fly.append(np.nan)
                pre_anti_saccade_peaks_per_fly.append(np.nan)
                continue
            else:
                pass

            syn_saccade_times, anti_saccade_times, anti_saccade_response_list, syn_saccade_response_list, pre_anti_saccade_response_list, pre_syn_saccade_response_list, \
                syn_saccade_turns, anti_saccade_turns, syn_saccade_peaks, anti_saccade_peaks, pre_syn_saccade_turns, pre_anti_saccade_turns, pre_syn_saccade_peaks, pre_anti_saccade_peaks\
                = time_between_saccades([[data[0][i][fly_indices]]], [saccade_turn_list[i][fly_indices]], p=p, start=0.0, duration=5.0)
            syn_saccade_response_list_per_fly.append(np.nanmean(syn_saccade_response_list))
            anti_saccade_response_list_per_fly.append(np.nanmean(anti_saccade_response_list))
            pre_anti_saccade_response_list_per_fly.append(np.nanmean(pre_anti_saccade_response_list))
            pre_syn_saccade_response_list_per_fly.append(np.nanmean(pre_syn_saccade_response_list))
            syn_saccade_turns_per_fly.append(np.nanmean(syn_saccade_turns))
            anti_saccade_turns_per_fly.append(np.nanmean(anti_saccade_turns))
            pre_syn_saccade_turns_per_fly.append(np.nanmean(pre_syn_saccade_turns))
            pre_anti_saccade_turns_per_fly.append(np.nanmean(pre_anti_saccade_turns))
            syn_saccade_peaks_per_fly.append(np.nanmean(syn_saccade_peaks))
            anti_saccade_peaks_per_fly.append(np.nanmean(anti_saccade_peaks))
            pre_syn_saccade_peaks_per_fly.append(np.nanmean(pre_syn_saccade_peaks))
            pre_anti_saccade_peaks_per_fly.append(np.nanmean(pre_anti_saccade_peaks))
        ##append values to the saccade_time_data_per_fly dict
        saccade_time_data_per_fly['syn_saccade_dist'].append(syn_saccade_response_list_per_fly)
        saccade_time_data_per_fly['anti_saccade_dist'].append(anti_saccade_response_list_per_fly)
        saccade_time_data_per_fly['syn_pre_saccade_dist'].append(pre_anti_saccade_response_list_per_fly)
        saccade_time_data_per_fly['anti_pre_saccade_dist'].append(pre_syn_saccade_response_list_per_fly)
        saccade_time_data_per_fly['syn_saccade_turns'].append(syn_saccade_turns_per_fly)
        saccade_time_data_per_fly['anti_saccade_turns'].append(anti_saccade_turns_per_fly)
        saccade_time_data_per_fly['syn_saccade_peaks'].append(syn_saccade_peaks_per_fly)
        saccade_time_data_per_fly['anti_saccade_peaks'].append(anti_saccade_peaks_per_fly)
        saccade_time_data_per_fly['pre_syn_saccade_turns'].append(syn_saccade_turns_per_fly)
        saccade_time_data_per_fly['pre_anti_saccade_turns'].append(anti_saccade_turns_per_fly)
        saccade_time_data_per_fly['pre_syn_saccade_peaks'].append(syn_saccade_peaks_per_fly)
        saccade_time_data_per_fly['pre_anti_saccade_peaks'].append(anti_saccade_peaks_per_fly)

        syn_saccade_turns_per_fly_all_stimuli.append(syn_saccade_turns_per_fly)
        anti_saccade_turns_per_fly_all_stimuli.append(anti_saccade_turns_per_fly)
        syn_saccade_peaks_per_fly_all_stimuli.append(syn_saccade_peaks_per_fly)
        anti_saccade_peaks_per_fly_all_stimuli.append(anti_saccade_peaks_per_fly)
        pre_syn_saccade_turns_per_fly_all_stimuli.append(pre_syn_saccade_turns_per_fly)
        pre_anti_saccade_turns_per_fly_all_stimuli.append(pre_anti_saccade_turns_per_fly)
        pre_syn_saccade_peaks_per_fly_all_stimuli.append(pre_syn_saccade_peaks_per_fly)
        pre_anti_saccade_peaks_per_fly_all_stimuli.append(pre_anti_saccade_peaks_per_fly)

    saccade_time_data_per_fly['syn_saccade_turns_allstimuli']=np.mean(syn_saccade_turns_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['anti_saccade_turns_allstimuli']=np.mean(anti_saccade_turns_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['syn_saccade_peaks_allstimuli']=np.mean(syn_saccade_peaks_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['anti_saccade_peaks_allstimuli']=np.mean(anti_saccade_peaks_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['pre_syn_saccade_turns_allstimuli']=np.mean(pre_syn_saccade_turns_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['pre_anti_saccade_turns_allstimuli']=np.mean(pre_anti_saccade_turns_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['pre_syn_saccade_peaks_allstimuli']=np.mean(pre_syn_saccade_peaks_per_fly_all_stimuli, axis=0)
    saccade_time_data_per_fly['pre_anti_saccade_peaks_allstimuli']=np.mean(pre_anti_saccade_peaks_per_fly_all_stimuli, axis=0)

    np.save(os.path.join(Dir_processed_data, filename + '_timewise_saccade_hemi_per_fly.npy'), saccade_time_data_per_fly, allow_pickle=True)

    color = [color_dict['FlpND'][1]]
    saccadic_and_smooth_response = {'variables':labels, 'conditions':conditions, 'syn_saccade_removed_response':[], 'syn_saccade_removed_angular_speed':[],
                                    'anti_saccade_removed_response': [], 'anti_saccade_removed_angular_speed': [],
                                    'all_saccade_removed_response': [], 'all_saccade_removed_angular_speed': [],
                                    'saccadic_response': [], 'saccadic_angular_speed': [], 'total_response': [], 'total_angular_speed':[]}
    opto_response_total = [opto_response_list[indices[i]] for i in range(len(indices))]
    angular_speed_total = [opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]
    saccadic_and_smooth_response['total_response'] = [np.mean(opto_response_list[indices[i]], axis=0) for i in range(len(indices))]
    saccadic_and_smooth_response['total_angular_speed'] = [np.mean(opto_response_angular_speed_list[indices[i]], axis=0) for i in range(len(indices))]


    # ## keep only anti-saccadic responses
    new_optomotor_response, new_opto_response_angular_speed = keep_saccades(opto_response_angular_speed, optomotor_response,
                            saccade_peak, window=[5, 5], p=p, keep='anti', stim_period = 0)
    ## multiply with expected response direction to fix saccade sign
    new_opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=new_optomotor_response[:, p+4:])
    new_opto_response_list = np.multiply(new_opto_response_list, saccade_peak[:, 0].reshape(-1, 1) * -1)
    new_opto_response_angular_speed_list = np.multiply(new_opto_response_angular_speed[:, p+3:], saccade_peak[:, 0].reshape(-1, 1) * -1)
    angular_speed_anti_saccades = [new_opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]
    opto_response_anti_saccades = [new_opto_response_list[indices[i]] for i in range(len(indices))]
    saccadic_and_smooth_response['anti_saccade_response'] = [np.mean(new_opto_response_list[indices[i]], axis=0) for i in range(len(indices))]
    saccadic_and_smooth_response['anti_saccade_angular_speed'] = [np.mean(new_opto_response_angular_speed_list[indices[i]], axis=0) for i in range(len(indices))]

    ## keep only syn-saccadic responses
    new_optomotor_response, new_opto_response_angular_speed = keep_saccades(opto_response_angular_speed, optomotor_response,
                            saccade_peak, window=[5, 5], p=p, keep='syn', stim_period = 0)
    ## multiply with expected response direction to fix saccade sign
    new_opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=new_optomotor_response[:, p+4:])
    new_opto_response_list = np.multiply(new_opto_response_list, saccade_peak[:, 0].reshape(-1, 1) * -1)
    new_opto_response_angular_speed_list = np.multiply(new_opto_response_angular_speed[:, p+3:], saccade_peak[:, 0].reshape(-1, 1) * -1)
    angular_speed_syn_saccades = [new_opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]
    opto_response_syn_saccades = [new_opto_response_list[indices[i]] for i in range(len(indices))]
    saccadic_and_smooth_response['syn_saccade_response'] = [np.mean(new_opto_response_list[indices[i]], axis=0) for i in range(len(indices))]
    saccadic_and_smooth_response['syn_saccade_angular_speed'] = [np.mean(new_opto_response_angular_speed_list[indices[i]], axis=0) for i in range(len(indices))]

    ## remove all saccades from the responses
    new_optomotor_response, new_opto_response_angular_speed = remove_saccades(opto_response_angular_speed, optomotor_response,
                            saccade_peak, window=[5, 5], p=p, remove='all', stim_period = 0, filler_type='zero')
    ## multiply with expected response direction to fix saccade sign
    new_opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=new_optomotor_response[:, p+4:])
    new_opto_response_list = np.multiply(new_opto_response_list, saccade_peak[:, 0].reshape(-1, 1) * -1)
    new_opto_response_angular_speed_list = np.multiply(new_opto_response_angular_speed[:, p+3:], saccade_peak[:, 0].reshape(-1, 1) * -1)
    opto_response_all_saccades_removed = [new_opto_response_list[indices[i]] for i in range(len(indices))]
    angular_speed_all_saccades_removed = [new_opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]
    saccadic_and_smooth_response['all_saccade_removed_response'] = [np.mean(new_opto_response_list[indices[i]], axis=0) for i in range(len(indices))]
    saccadic_and_smooth_response['all_saccade_removed_angular_speed'] = [np.mean(new_opto_response_angular_speed_list[indices[i]], axis=0) for i in range(len(indices))]

    ## remove all saccades from the responses, linearly
    new_optomotor_response, new_opto_response_angular_speed = remove_saccades(opto_response_angular_speed, optomotor_response,
                            saccade_peak, window=[5, 5], p=p, remove='all', stim_period = 0, filler_type='linear')
    ## multiply with expected response direction to fix saccade sign
    new_opto_response_angular_speed_list = np.multiply(new_opto_response_angular_speed[:, p+3:], saccade_peak[:, 0].reshape(-1, 1) * -1)
    angular_speed_all_saccades_removed_linear = [new_opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]
    plot_power_spectrum([opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]+angular_speed_all_saccades_removed_linear,
                        Dir, filename, ['r', 'g', 'grey', 'r', 'g', 'grey'], ls=['solid', 'solid', 'solid', 'dashed', 'dashed', 'dashed'])

    ## keep only saccadic responses
    new_optomotor_response, new_opto_response_angular_speed = keep_saccades(opto_response_angular_speed, optomotor_response, saccade_peak, window=[5, 5], p=p, keep='all', stim_period = 0)
    ## multiply with expected response direction to fix saccade sign
    new_opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=new_optomotor_response[:, p+4:])
    new_opto_response_list = np.multiply(new_opto_response_list, saccade_peak[:, 0].reshape(-1, 1) * -1)
    new_opto_response_angular_speed_list = np.multiply(new_opto_response_angular_speed[:, p+3:], saccade_peak[:, 0].reshape(-1, 1) * -1)
    opto_response_only_saccades = [new_opto_response_list[indices[i]] for i in range(len(indices))]
    angular_speed_only_saccades = [new_opto_response_angular_speed_list[indices[i]] for i in range(len(indices))]
    saccadic_and_smooth_response['all_saccade_response'] = [np.mean(new_opto_response_list[indices[i]], axis=0) for i in range(len(indices))]
    saccadic_and_smooth_response['all_saccade_angular_speed'] = [np.mean(new_opto_response_angular_speed_list[indices[i]], axis=0) for i in range(len(indices))]


    ## prepare data for each fly
    fly_id = optomotor_response[:, p + 2]
    fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
    saccadic_angular_speed_list = []
    smooth_angular_speed_list = []
    saccadic_response_list = []
    smooth_response_list = []
    syn_saccadic_response_list = []
    anti_saccadic_response_list = []

    saccadic_response_eachfly_list = []
    smooth_response_eachfly_list = []
    saccadic_angular_speed_eachfly_list = []
    smooth_angular_speed_eachfly_list = []
    total_angular_speed_eachfly_list = []
    syn_saccadic_response_eachfly_list = []
    anti_saccadic_response_eachfly_list = []

    saccadic_angular_speed_sem_eachfly_list = []
    smooth_angular_speed_sem_eachfly_list = []
    syn_saccadic_response_sem_eachfly_list = []
    anti_saccadic_response_sem_eachfly_list = []
    for i in range(len(angular_speed_only_saccades)):
        saccadic_angular_speed_eachfly = []
        smooth_angular_speed_eachfly = []
        total_angular_speed_eachfly = []
        saccadic_response_eachfly = []
        smooth_response_eachfly = []
        syn_saccadic_response_eachfly = []
        anti_saccadic_response_eachfly = []

        saccadic_angular_speed_sem_eachfly = []
        smooth_angular_speed_sem_eachfly = []
        syn_saccadic_response_sem_eachfly = []
        anti_saccadic_response_sem_eachfly = []

        # all
        saccadic_angular_speed_list.append(np.mean(angular_speed_only_saccades[i][:, stim_period:], axis=1))
        smooth_angular_speed_list.append(np.mean(angular_speed_all_saccades_removed[i][:, stim_period:], axis=1))
        saccadic_response_list.append(opto_response_only_saccades[i][:, -1])
        smooth_response_list.append(opto_response_all_saccades_removed[i][:, -1])
        syn_saccadic_response_list.append(opto_response_syn_saccades[i][:, -1])
        anti_saccadic_response_list.append(opto_response_anti_saccades[i][:, -1])

        for ids in np.unique(fly_id_indices[i]):
            fly_indices = np.argwhere(fly_id_indices[i] == ids).flatten()
            if fly_indices.shape[0] < 5: #this value is either 10 or 5
                continue
            else:
                # mean
                total_angular_speed_eachfly.append(np.mean(angular_speed_total[i][fly_indices, stim_period:]))
                saccadic_angular_speed_eachfly.append(np.mean(angular_speed_only_saccades[i][fly_indices, stim_period:]))
                smooth_angular_speed_eachfly.append(np.mean(angular_speed_all_saccades_removed[i][fly_indices, stim_period:]))
                saccadic_response_eachfly.append(np.mean(opto_response_only_saccades[i][fly_indices, stim_period:]))
                smooth_response_eachfly.append(np.mean(opto_response_all_saccades_removed[i][fly_indices, stim_period:]))
                syn_saccadic_response_eachfly.append(np.mean(opto_response_syn_saccades[i][fly_indices, -1]))
                anti_saccadic_response_eachfly.append(np.mean(opto_response_anti_saccades[i][fly_indices, -1]))
                # sem
                saccadic_angular_speed_sem_eachfly.append(sem(np.mean(angular_speed_only_saccades[i][fly_indices, stim_period:], axis=1)))
                smooth_angular_speed_sem_eachfly.append(sem(np.mean(angular_speed_all_saccades_removed[i][fly_indices, stim_period:], axis=1)))
                syn_saccadic_response_sem_eachfly.append(sem(opto_response_syn_saccades[i][fly_indices, -1]))
                anti_saccadic_response_sem_eachfly.append(sem(opto_response_anti_saccades[i][fly_indices, -1]))
        saccadic_angular_speed_eachfly_list.append(saccadic_angular_speed_eachfly)
        smooth_angular_speed_eachfly_list.append(smooth_angular_speed_eachfly)
        saccadic_response_eachfly_list.append(saccadic_response_eachfly)
        smooth_response_eachfly_list.append(smooth_response_eachfly)
        syn_saccadic_response_eachfly_list.append(syn_saccadic_response_eachfly)
        anti_saccadic_response_eachfly_list.append(anti_saccadic_response_eachfly)
        total_angular_speed_eachfly_list.append(total_angular_speed_eachfly)

        saccadic_angular_speed_sem_eachfly_list.append(saccadic_angular_speed_sem_eachfly)
        smooth_angular_speed_sem_eachfly_list.append(smooth_angular_speed_sem_eachfly)
        syn_saccadic_response_sem_eachfly_list.append(syn_saccadic_response_sem_eachfly)
        anti_saccadic_response_sem_eachfly_list.append(anti_saccadic_response_sem_eachfly)

    saccadic_and_smooth_response['saccadic_response_all'] = saccadic_angular_speed_list
    saccadic_and_smooth_response['smooth_response_all'] = smooth_angular_speed_list
    saccadic_and_smooth_response['saccadic_turns_all'] = saccadic_response_list
    saccadic_and_smooth_response['smooth_turns_all'] = smooth_response_list
    saccadic_and_smooth_response['syn_saccadic_turns_all'] = syn_saccadic_response_list
    saccadic_and_smooth_response['anti_saccadic_turns_all'] = anti_saccadic_response_list

    saccadic_and_smooth_response['saccadic_response_eachfly'] = saccadic_angular_speed_eachfly_list
    saccadic_and_smooth_response['smooth_response_eachfly'] = smooth_angular_speed_eachfly_list
    saccadic_and_smooth_response['saccadic_turns_eachfly'] = saccadic_response_eachfly_list
    saccadic_and_smooth_response['smooth_turns_eachfly'] = smooth_response_eachfly_list
    saccadic_and_smooth_response['syn_saccadic_turns_eachfly'] = syn_saccadic_response_eachfly_list
    saccadic_and_smooth_response['anti_saccadic_turns_eachfly'] = anti_saccadic_response_eachfly_list
    saccadic_and_smooth_response['total_response_eachfly'] = total_angular_speed_eachfly_list

    saccadic_and_smooth_response['saccadic_response_sem_eachfly'] = saccadic_angular_speed_sem_eachfly_list
    saccadic_and_smooth_response['smooth_response_sem_eachfly'] = smooth_angular_speed_sem_eachfly_list
    saccadic_and_smooth_response['syn_saccadic_sem_turns_eachfly'] = syn_saccadic_response_sem_eachfly_list
    saccadic_and_smooth_response['anti_saccadic_sem_turns_eachfly'] = anti_saccadic_response_sem_eachfly_list

    ## plot only saccadic response
    saccadic_and_smooth_response['saccadic_response'] = [np.mean(opto_response_list[indices[i]], axis=0) - np.mean(new_opto_response_list[indices[i]], axis=0) for i in range(len(indices))]
    saccadic_and_smooth_response['saccadic_angular_speed'] = [np.mean(opto_response_angular_speed_list[indices[i]], axis=0) - np.mean(new_opto_response_angular_speed_list[indices[i]], axis=0)
                                                              for i in range(len(indices))]
    np.save(os.path.join(Dir_processed_data, filename+'_saccadic_smooth_response.npy'), saccadic_and_smooth_response, allow_pickle=True)
    ## linear comparison for saccadic and smooth response separately
    ## saccadic
    _, _, _, _, _, saccadic_cumulative_error_list, saccadic_prediction_error = draw_trace_plot(angular_speed_only_saccades, Dir, filename + 'saccadic_angular_speed_linear_comparison',
                labels=labels, top=150, bottom=-50, cmap=ListedColormap(color * 3), smoothed='Y', highlight=[0, 0], ls=['dotted', 'dashed', 'solid'],
                draw_trace='N', draw_linear_sum='Y', draw_median='N', start=180, end=-1, axvline=120, fly_indices=fly_id_indices, PDF='Y')
    ## smooth
    _, _, _, _, _, smooth_cumulative_error_list, smooth_prediction_error = draw_trace_plot(angular_speed_all_saccades_removed, Dir, filename + 'smooth_angular_speed_linear_comparison',
                labels=labels, top=150, bottom=-50, cmap=ListedColormap(color * 3), smoothed='Y', highlight=[0, 0], ls=['dotted', 'dashed', 'solid'],
                draw_trace='N', draw_linear_sum='Y', draw_median='N', start=180, end=-1, axvline=120, fly_indices=fly_id_indices, PDF='Y')

    total_data = {'cumulative_error_list': smooth_cumulative_error_list, 'prediction_error': smooth_prediction_error}
    np.save(os.path.join(Dir_processed_data, filename + '_smooth_cumulative_error_list.npy'), total_data, allow_pickle=True)
    total_data = {'cumulative_error_list': saccadic_cumulative_error_list, 'prediction_error': saccadic_prediction_error}
    np.save(os.path.join(Dir_processed_data, filename + '_saccadic_cumulative_error_list.npy'), total_data, allow_pickle=True)

    cw_anti, ccw_anti, cw_syn, ccw_syn = extract_saccades(optomotor_response, opto_response_angular_speed, saccade_peak, indices, window=[10, 10], p=p)
    for i in range(len(cw_syn)):
        draw_trace_plot([cw_syn[i][0], ccw_anti[i][0], ccw_syn[i][0], cw_anti[i][0]], Dir, filename+'saccade_comparison', ['cw_syn', 'ccw_anti', 'ccw_syn', 'cw_anti'], top=60,
                        bottom=-60, cmap='Accent', smoothed='N', highlight=None, filter='N')
        draw_trace_plot([cw_syn[i][1], ccw_anti[i][1], ccw_syn[i][1], cw_anti[i][1]], Dir, filename + 'saccade_comparison_angular_speed', ['cw_syn', 'ccw_anti', 'ccw_syn', 'cw_anti'], top=1000,
                        bottom=-1000, cmap='Accent', smoothed='N', highlight=None, filter='N')
    saccade_shapes = {'variables':labels, 'conditions':conditions, 'cw_syn_angular_speed': [np.mean(cw_syn[i][0], axis=0) for i in range(len(indices))],
                      'ccw_syn_angular_speed': [np.mean(ccw_syn[i][0], axis=0) for i in range(len(indices))], 'cw_anti_angular_speed': [np.mean(cw_anti[i][0], axis=0) for i in range(len(indices))],
                      'ccw_anti_angular_speed': [np.mean(ccw_anti[i][0], axis=0) for i in range(len(indices))], 'cw_syn_opto_response': [np.mean(cw_syn[i][1], axis=0) for i in range(len(indices))],
                      'ccw_syn_opto_response': [np.mean(ccw_syn[i][1], axis=0) for i in range(len(indices))], 'cw_anti_opto_response': [np.mean(cw_anti[i][1], axis=0) for i in range(len(indices))],
                      'ccw_anti_opto_response': [np.mean(ccw_anti[i][1], axis=0) for i in range(len(indices))]}
    np.save(os.path.join(Dir_processed_data, filename+'_saccade_shape_hemi.npy'), saccade_shapes, allow_pickle=True)
    return 1


def syn_anti_saccade_distribution(ccw_syn_anti_saccades_per_trial, cw_syn_anti_saccades_per_trial, Dir, filename):
    fig, ax = plt.subplots()
    labels = [2, 5, 10, 15, 20, 25]
    fig, ax = plt.subplots(2, 3)
    m = 3

    for i in range(len(ccw_syn_anti_saccades_per_trial)):
        heatmap, xedges, yedges = np.histogram2d(ccw_syn_anti_saccades_per_trial[i][0, :], ccw_syn_anti_saccades_per_trial[i][1, :], bins=(29, 29), range=[[0, 30], [0, 30]],
                                                 density=True)
        im = ax[i // m, i % m].imshow(heatmap, origin='lower', extent=[0, 30, 0, 30], cmap='afmhot', vmin=0, vmax=0.1)
        for pos in ['left', 'right', 'top', 'bottom']:
            ax[i // m, i % m].spines[pos].set_visible(False)
        ax[i // m, i % m].set_xticks([])
        ax[i // m, i % m].set_yticks([])
    for i in range(len(cw_syn_anti_saccades_per_trial)):
        heatmap, xedges, yedges = np.histogram2d(cw_syn_anti_saccades_per_trial[i][0, :], cw_syn_anti_saccades_per_trial[i][1, :], bins=(29, 29), range=[[0, 30], [0, 30]],
                                                 density=True)
        im = ax[i // m, i % m].imshow(heatmap, origin='lower', extent=[0, 30, 0, 30], cmap='afmhot', vmin=0, vmax=0.1)
        for pos in ['left', 'right', 'top', 'bottom']:
            ax[i // m, i % m].spines[pos].set_visible(False)
        ax[i // m, i % m].set_xticks([])
        ax[i // m, i % m].set_yticks([])
    fig.savefig(os.path.join(Dir, filename+'syn_anti_dist.png'))
    return 1


def draw_saccade_size_comparison(Dir, filename, data, labels, cmap='tab10'):
    ## data contains opto_response, angular_speed, saccade_peaks divide by conditions
    cw_anti_list = []
    ccw_anti_list = []
    cw_syn_list = []
    ccw_syn_list = []
    for i in range(len(data[0])):
        cw_anti, ccw_anti, cw_syn, ccw_syn = extract_saccades(data[0][i], data[1][i], data[2][i], window=[3, 5])
        ccw_anti_list.append(ccw_anti[0])
        cw_anti_list.append(cw_anti[0])
        cw_syn_list.append(cw_syn[0])
        ccw_syn_list.append(ccw_syn[0])

    n = len(data[0])
    if cm.get_cmap(cmap).N == 256:
        colors = [cm.get_cmap(cmap)(i) for i in np.linspace(0, 1, n)]
    else:
        colors = [cm.get_cmap(cmap)(i) for i in np.arange(0, n)]

    fig, ax = plt.subplots(1, 4, figsize=(figsize[0]*4, figsize[1]))
    for i in range(len(cw_anti_list)):
        ax[0].plot(np.mean(cw_anti_list[i], axis=0), c=colors[i], lw=3, label=labels[i])
        ax[1].plot(np.mean(ccw_anti_list[i], axis=0), c=colors[i], lw=3, label=labels[i])
        ax[2].plot(np.mean(cw_syn_list[i], axis=0), c=colors[i], lw=3, label=labels[i])
        ax[3].plot(np.mean(ccw_syn_list[i], axis=0), c=colors[i], lw=3, label=labels[i])
    ax[3].legend(loc='upper right')
    ax[0].set_title('cw_anti')
    ax[1].set_title('ccw_anti')
    ax[2].set_title('cw_syn')
    ax[3].set_title('ccw_syn')
    ax[1].set_yticklabels([])
    ax[2].set_yticklabels([])
    ax[3].set_yticklabels([])
    for i in range(4):
        ax[i].set_ylim(top=80, bottom=-80)
    fig.savefig(os.path.join(Dir, filename+'_traces.png'))

    cw_anti_size_list = []
    ccw_anti_size_list = []
    cw_syn_size_list = []
    ccw_syn_size_list = []
    for i in range(len(data[0])):
        ccw_anti_size_list.append(np.apply_along_axis(lambda x:np.subtract(x[-1], x[0]), arr=ccw_anti_list[i], axis= 1))
        cw_anti_size_list.append(np.apply_along_axis(lambda x:np.subtract(x[-1], x[0]), arr=cw_anti_list[i], axis= 1))
        cw_syn_size_list.append(np.apply_along_axis(lambda x:np.subtract(x[-1], x[0]), arr=cw_syn_list[i], axis= 1))
        ccw_syn_size_list.append(np.apply_along_axis(lambda x:np.subtract(x[-1], x[0]), arr=ccw_syn_list[i], axis= 1))
    fig, ax = plt.subplots(1, 4, figsize=(figsize[0]*4, figsize[1]))
    positions = np.linspace(0,n-1,n)
    result = ax[0].violinplot(cw_anti_size_list, positions, showmeans=True, showextrema=False, widths=0.5)
    for i in range(len(result['bodies'])):
        result['bodies'][i].set_facecolor((0, 1, 0, 0.6))
    result = ax[1].violinplot(ccw_anti_size_list, positions, showmeans=True, showextrema=False, widths=0.5)
    for i in range(len(result['bodies'])):
        result['bodies'][i].set_facecolor((0, 1, 0, 0.6))
    result = ax[2].violinplot(cw_syn_size_list, positions, showmeans=True, showextrema=False, widths=0.5)
    for i in range(len(result['bodies'])):
        result['bodies'][i].set_facecolor((0, 1, 0, 0.6))
    result = ax[3].violinplot(ccw_syn_size_list, positions, showmeans=True, showextrema=False, widths=0.5)
    for i in range(len(result['bodies'])):
        result['bodies'][i].set_facecolor((0, 1, 0, 0.6))
    ax[0].set_title('cw_anti')
    ax[1].set_title('ccw_anti')
    ax[2].set_title('cw_syn')
    ax[3].set_title('ccw_syn')
    ax[0].set_ylim(top=100, bottom=-100)
    ax[1].set_ylim(top=100, bottom=-100)
    ax[2].set_ylim(top=100, bottom=-100)
    ax[3].set_ylim(top=100, bottom=-100)
    ax[1].set_yticklabels([])
    ax[2].set_yticklabels([])
    ax[3].set_yticklabels([])
    fig.savefig(os.path.join(Dir, filename + '_violin.png'))
    return cw_anti_list


def draw_saccade_data(data, Dir, filename, labels, p=2):
    syn_saccade_responses_list = []
    anti_saccade_responses_list = []
    pre_syn_saccade_responses_list = []
    pre_anti_saccade_responses_list = []
    trials = []
    for i in range(len(data[0])):
        cw_indices = np.argwhere(data[0][i][:, 0] == -1).flatten()
        cw_syn_indices = np.argwhere(data[0][i][cw_indices][:, p+2+stim_period:] > 0)
        cw_anti_indices = np.argwhere(data[0][i][cw_indices][:, p+2+stim_period:] < 0)

        ccw_indices = np.argwhere(data[0][i][:, 0] == 1).flatten()
        ccw_syn_indices = np.argwhere(data[0][i][ccw_indices][:, p+2+stim_period:] < 0)
        ccw_anti_indices = np.argwhere(data[0][i][ccw_indices][:, p+2+stim_period:] > 0)

        syn_saccade_responses_list.append(cw_syn_indices.shape[0] + ccw_syn_indices.shape[0])
        anti_saccade_responses_list.append(cw_anti_indices.shape[0] + ccw_anti_indices.shape[0])

        cw_indices = np.argwhere(data[0][i][:, 0] == -1).flatten()
        cw_syn_indices = np.argwhere(data[0][i][cw_indices][:, p+2:stim_period+p+2] > 0)
        cw_anti_indices = np.argwhere(data[0][i][cw_indices][:, p+2:stim_period+p+2] < 0)

        ccw_indices = np.argwhere(data[0][i][:, 0] == 1).flatten()
        ccw_syn_indices = np.argwhere(data[0][i][ccw_indices][:, p+2:stim_period+p+2] < 0)
        ccw_anti_indices = np.argwhere(data[0][i][ccw_indices][:, p+2:stim_period+p+2] > 0)

        pre_syn_saccade_responses_list.append(cw_syn_indices.shape[0] + ccw_syn_indices.shape[0])
        pre_anti_saccade_responses_list.append(cw_anti_indices.shape[0] + ccw_anti_indices.shape[0])
        trials.append(data[0][i].shape[0])

    draw_bar_plot_saccade(trials, anti_saccade_responses_list, syn_saccade_responses_list, Dir, filename+'_200_saccade_responsenew', labels)
    draw_bar_plot_saccade(trials, pre_anti_saccade_responses_list, pre_syn_saccade_responses_list, Dir, filename+'_pre_saccade_responsenew', labels)
    draw_bar_plot_saccade(trials, np.subtract(anti_saccade_responses_list, pre_anti_saccade_responses_list), np.subtract(syn_saccade_responses_list,
                                            pre_syn_saccade_responses_list),Dir, filename+'_change_in_saccades', labels, ylim=[10, -2])
    return 1


def draw_bar_plot_saccade(anti_saccade_responses_list, syn_saccade_responses_list, Dir, filename, labels, ylim=[15, 0], color=[]):
    # cmap = ListedColormap(all_color_data['RdBu3'])
    # colors = [cmap(0.9), cmap(0.1)]
    colors = [(122/255, 204/255, 230/255), (9/255, 73/255, 145/255)]
    # colors = [(230 / 255, 181 / 255, 210 / 255), (191 / 255, 85 / 255, 105 / 255)]
    draw_bar_plot([syn_saccade_responses_list, anti_saccade_responses_list], Dir, filename, labels, ylim=ylim, colors=colors, names=['Anti', 'Syn'], PDF='Y')
    # draw_violin_plot([syn_saccade_responses_list, anti_saccade_responses_list], ['Anti', 'Syn'], labels, Dir, filename, ylim=ylim, colors=colors)
    return 1


def time_between_saccades(data, saccade_turns, p=2, start=0.0, duration=2.5):
    # data, conditions, indices = split_data_by_conditions([saccade_peak], [saccade_peak[:, 1]])

    start = int(start*60)
    duration = int(duration*60)
    syn_saccade_times = []
    anti_saccade_times = []
    cw_syn_anti_saccades_per_trial = []
    ccw_syn_anti_saccades_per_trial = []
    syn_saccade_turns = []
    anti_saccade_turns = []
    syn_saccade_peaks = []
    anti_saccade_peaks = []
    for i in range(len(data[0])):
        cw_indices = np.argwhere(data[0][i][:, 0] == -1).flatten()
        cw_syn_indices = np.argwhere(data[0][i][cw_indices][:, p+2+stim_period+start:p+2+stim_period+start+duration] > 0)
        cw_syn_saccades = np.diff(cw_syn_indices[:, 1])
        cw_syn_saccades = cw_syn_saccades[cw_syn_saccades > 0]
        no_of_cw_syn_saccades_per_trial = np.array([np.argwhere(cw_syn_indices[:, 0]==x).flatten().shape[0] for x in range(cw_indices.shape[0])])
        cw_anti_indices = np.argwhere(data[0][i][cw_indices][:, p+2+stim_period+start:p+2+stim_period+start+duration] < 0)
        cw_anti_saccades = np.diff(cw_anti_indices[:, 1])
        cw_anti_saccades = cw_anti_saccades[cw_anti_saccades > 0]
        no_of_cw_anti_saccades_per_trial = np.array([np.argwhere(cw_anti_indices[:, 0] == x).flatten().shape[0] for x in range(cw_indices.shape[0])])

        ccw_indices = np.argwhere(data[0][i][:, 0] == 1).flatten()
        ccw_syn_indices = np.argwhere(data[0][i][ccw_indices][:, p+2+stim_period+start:p+2+stim_period+start+duration] < 0)
        ccw_syn_saccades = np.diff(ccw_syn_indices[:, 1])
        ccw_syn_saccades = ccw_syn_saccades[ccw_syn_saccades > 0]
        no_of_ccw_syn_saccades_per_trial = np.array([np.argwhere(ccw_syn_indices[:, 0] == x).flatten().shape[0] for x in range(ccw_indices.shape[0])])
        ccw_anti_indices = np.argwhere(data[0][i][ccw_indices][:, p+2+stim_period+start:p+2+stim_period+start+duration] > 0)
        ccw_anti_saccades = np.diff(ccw_anti_indices[:, 1])
        ccw_anti_saccades = ccw_anti_saccades[ccw_anti_saccades > 0]
        no_of_ccw_anti_saccades_per_trial = np.array([np.argwhere(ccw_anti_indices[:, 0] == x).flatten().shape[0] for x in range(ccw_indices.shape[0])])

        syn_saccade_turns.append(np.concatenate((saccade_turns[i][ccw_indices][(ccw_syn_indices[:, 0], p + 2 + stim_period + start + ccw_syn_indices[:, 1])],
                                                saccade_turns[i][cw_indices][(cw_syn_indices[:, 0], p+2+stim_period+start+cw_syn_indices[:, 1])]), axis=None))
        anti_saccade_turns.append(np.concatenate((saccade_turns[i][ccw_indices][(ccw_anti_indices[:, 0], p + 2 + stim_period + start + ccw_anti_indices[:, 1])],
                                                saccade_turns[i][cw_indices][(cw_anti_indices[:, 0], p+2+stim_period+start+cw_anti_indices[:, 1])]), axis=None))

        syn_saccade_peaks.append(np.concatenate((data[0][i][ccw_indices][(ccw_syn_indices[:, 0], p + 2 + stim_period + start + ccw_syn_indices[:, 1])]*-1,
                                                data[0][i][cw_indices][(cw_syn_indices[:, 0], p+2+stim_period+start+cw_syn_indices[:, 1])]), axis=None))
        anti_saccade_peaks.append(np.concatenate((data[0][i][ccw_indices][(ccw_anti_indices[:, 0], p + 2 + stim_period + start + ccw_anti_indices[:, 1])]*-1,
                                                data[0][i][cw_indices][(cw_anti_indices[:, 0], p+2+stim_period+start+cw_anti_indices[:, 1])]), axis=None))

        syn_saccade_times.append(np.concatenate((cw_syn_saccades, ccw_syn_saccades)))
        anti_saccade_times.append(np.concatenate((cw_anti_saccades, ccw_anti_saccades)))

        cw_syn_anti_saccades_per_trial.append(np.array([no_of_cw_syn_saccades_per_trial, no_of_cw_anti_saccades_per_trial]))
        ccw_syn_anti_saccades_per_trial.append(np.array([no_of_ccw_syn_saccades_per_trial, no_of_ccw_anti_saccades_per_trial]))

    anti_saccade_response_list = []
    syn_saccade_response_list = []
    for i in range(len(cw_syn_anti_saccades_per_trial)):
        anti_saccade_response_list.append(np.concatenate((cw_syn_anti_saccades_per_trial[i][1], ccw_syn_anti_saccades_per_trial[i][1])))
        syn_saccade_response_list.append(np.concatenate((cw_syn_anti_saccades_per_trial[i][0], ccw_syn_anti_saccades_per_trial[i][0])))

    pre_syn_saccade_times = []
    pre_anti_saccade_times = []
    pre_cw_syn_anti_saccades_per_trial = []
    pre_ccw_syn_anti_saccades_per_trial = []
    pre_syn_saccade_turns = []
    pre_anti_saccade_turns = []
    pre_syn_saccade_peaks = []
    pre_anti_saccade_peaks = []
    for i in range(len(data[0])):
        cw_indices = np.argwhere(data[0][i][:, 0] == -1).flatten()
        cw_syn_indices = np.argwhere(data[0][i][cw_indices][:, p+2+(stim_period-duration):p+2+stim_period] > 0)
        cw_syn_saccades = np.diff(cw_syn_indices[:, 1])
        cw_syn_saccades = cw_syn_saccades[cw_syn_saccades > 0]
        no_of_cw_syn_saccades_per_trial = np.array([np.argwhere(cw_syn_indices[:, 0]==x).flatten().shape[0] for x in range(cw_indices.shape[0])])
        cw_anti_indices = np.argwhere(data[0][i][cw_indices][:, p+2+(stim_period-duration):p+2+stim_period] < 0)
        cw_anti_saccades = np.diff(cw_anti_indices[:, 1])
        cw_anti_saccades = cw_anti_saccades[cw_anti_saccades > 0]
        no_of_cw_anti_saccades_per_trial = np.array([np.argwhere(cw_anti_indices[:, 0] == x).flatten().shape[0] for x in range(cw_indices.shape[0])])

        ccw_indices = np.argwhere(data[0][i][:, 0] == 1).flatten()
        ccw_syn_indices = np.argwhere(data[0][i][ccw_indices][:, p+2+(stim_period-duration):p+2+stim_period] < 0)
        ccw_syn_saccades = np.diff(ccw_syn_indices[:, 1])
        ccw_syn_saccades = ccw_syn_saccades[ccw_syn_saccades > 0]
        no_of_ccw_syn_saccades_per_trial = np.array([np.argwhere(ccw_syn_indices[:, 0] == x).flatten().shape[0] for x in range(ccw_indices.shape[0])])
        ccw_anti_indices = np.argwhere(data[0][i][ccw_indices][:, p+2+(stim_period-duration):p+2+stim_period] > 0)
        ccw_anti_saccades = np.diff(ccw_anti_indices[:, 1])
        ccw_anti_saccades = ccw_anti_saccades[ccw_anti_saccades > 0]
        no_of_ccw_anti_saccades_per_trial = np.array([np.argwhere(ccw_anti_indices[:, 0] == x).flatten().shape[0] for x in range(ccw_indices.shape[0])])

        pre_syn_saccade_times.append(np.concatenate((cw_syn_saccades, ccw_syn_saccades)))
        pre_anti_saccade_times.append(np.concatenate((cw_anti_saccades, ccw_anti_saccades)))

        pre_cw_syn_anti_saccades_per_trial.append(np.array([no_of_cw_syn_saccades_per_trial, no_of_cw_anti_saccades_per_trial]))
        pre_ccw_syn_anti_saccades_per_trial.append(np.array([no_of_ccw_syn_saccades_per_trial, no_of_ccw_anti_saccades_per_trial]))

        pre_syn_saccade_turns.append(np.concatenate((saccade_turns[i][ccw_indices][(ccw_syn_indices[:, 0], p+2+(stim_period-duration) + ccw_syn_indices[:, 1])],
                                                 saccade_turns[i][cw_indices][(cw_syn_indices[:, 0], p+2+(stim_period-duration) + cw_syn_indices[:, 1])]), axis=None))
        pre_anti_saccade_turns.append(np.concatenate((saccade_turns[i][ccw_indices][(ccw_anti_indices[:, 0], p+2+(stim_period-duration) + ccw_anti_indices[:, 1])],
                                                  saccade_turns[i][cw_indices][(cw_anti_indices[:, 0], p+2+(stim_period-duration) + cw_anti_indices[:, 1])]), axis=None))

        pre_syn_saccade_peaks.append(np.concatenate((data[0][i][ccw_indices][(ccw_syn_indices[:, 0], p+2+(stim_period-duration) + ccw_syn_indices[:, 1])]*-1,
                                                 data[0][i][cw_indices][(cw_syn_indices[:, 0], p+2+(stim_period-duration) + cw_syn_indices[:, 1])]), axis=None))
        pre_anti_saccade_peaks.append(np.concatenate((data[0][i][ccw_indices][(ccw_anti_indices[:, 0], p+2+(stim_period-duration) + ccw_anti_indices[:, 1])]*-1,
                                                  data[0][i][cw_indices][(cw_anti_indices[:, 0], p+2+(stim_period-duration) + cw_anti_indices[:, 1])]), axis=None))

    pre_anti_saccade_response_list = []
    pre_syn_saccade_response_list = []
    for i in range(len(cw_syn_anti_saccades_per_trial)):
        pre_anti_saccade_response_list.append(np.concatenate((pre_cw_syn_anti_saccades_per_trial[i][1], pre_ccw_syn_anti_saccades_per_trial[i][1])))
        pre_syn_saccade_response_list.append(np.concatenate((pre_cw_syn_anti_saccades_per_trial[i][0], pre_ccw_syn_anti_saccades_per_trial[i][0])))
    return syn_saccade_times, anti_saccade_times, anti_saccade_response_list, syn_saccade_response_list, pre_anti_saccade_response_list, pre_syn_saccade_response_list,\
            syn_saccade_turns, anti_saccade_turns, syn_saccade_peaks, anti_saccade_peaks,pre_syn_saccade_turns, pre_anti_saccade_turns, pre_syn_saccade_peaks, pre_anti_saccade_peaks


def extract_saccades(opto_response, opto_angular_speed, saccade_peaks, indices, window=[7,7], p=2):
    opto_response = opto_response[:, p+4:]
    opto_angular_speed = opto_angular_speed[:, p+3:]

    cw_anti = []
    cw_syn = []
    ccw_anti = []
    ccw_syn = []
    for k in indices:
        cw_anti_saccades = []
        ccw_anti_saccades = []
        cw_syn_saccades = []
        ccw_syn_saccades = []
        cw_anti_saccades_angular_speed = []
        ccw_anti_saccades_angular_speed = []
        cw_syn_saccades_angular_speed = []
        ccw_syn_saccades_angular_speed = []
        for i in k:
            if saccade_peaks[i, 0] == -1:
                anti_indices = np.argwhere(saccade_peaks[i, p+2 + stim_period:] < 0)
                syn_indices = np.argwhere(saccade_peaks[i, p+2 + stim_period:] > 0)
                for index in anti_indices:
                    if index[0] > window[0] and index[0] < stim_period - window[1]:
                        cw_anti_saccades.append(np.subtract(opto_response[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]], opto_response[i, index[0] + stim_period]))
                        cw_anti_saccades_angular_speed.append(opto_angular_speed[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]])
                for index in syn_indices:
                    if index[0] > window[0] and index[0] < stim_period - window[1]:
                        cw_syn_saccades.append(np.subtract(opto_response[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]], opto_response[i, index[0] + stim_period]))
                        cw_syn_saccades_angular_speed.append(opto_angular_speed[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]])
            elif saccade_peaks[i, 0] == 1:
                anti_indices = np.argwhere(saccade_peaks[i, p+2 + stim_period:] > 0)
                syn_indices = np.argwhere(saccade_peaks[i, p+2 + stim_period:] < 0)
                for index in anti_indices:
                    if index[0] > window[0] and index[0] < stim_period - window[1]:
                        ccw_anti_saccades.append(np.subtract(opto_response[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]], opto_response[i, index[0] + stim_period]))
                        ccw_anti_saccades_angular_speed.append(opto_angular_speed[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]])
                for index in syn_indices:
                    if index[0] > window[0] and index[0] < stim_period - window[1]:
                        ccw_syn_saccades.append(np.subtract(opto_response[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]], opto_response[i, index[0] + stim_period]))
                        ccw_syn_saccades_angular_speed.append(opto_angular_speed[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]])
        cw_anti.append(np.array([np.array(cw_anti_saccades), np.array(cw_anti_saccades_angular_speed)]))
        cw_syn.append(np.array([np.array(cw_syn_saccades), np.array(cw_syn_saccades_angular_speed)]))
        ccw_anti.append(np.array([np.array(ccw_anti_saccades), np.array(ccw_anti_saccades_angular_speed)]))
        ccw_syn.append(np.array([np.array(ccw_syn_saccades), np.array(ccw_syn_saccades_angular_speed)]))
    return cw_anti, ccw_anti, cw_syn, ccw_syn
