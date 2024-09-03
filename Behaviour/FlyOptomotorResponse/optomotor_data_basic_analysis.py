import os
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

# add all your local imports here
from reject_nonresponsive_trials import *
from split_optomotor_data_by_stimulus import split_data_by_conditions, indices_by_condition
import draw_trace_each_fly
import draw_optomotor_response
import get_expected_response_hemi


def analysis(Dir, filename=None, sex='', p=3, condn_indices=[2, 1, 0], ignore_values=[], variables=[]):
    if filename==None:
        filename = Dir.split(r'\'')[0].split('\\')[-1]
    for files in os.listdir(Dir):
        if files.endswith(filename+'_opto_response_'+sex+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response = np.load(path, allow_pickle=True)
        elif files.endswith(filename+'_opto_response_'+sex+'_speed'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response_speed = np.load(path,allow_pickle=True)
        elif files.endswith(filename+'_opto_response_'+sex+'_angular_speed'+'.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response_angular_speed = np.load(path, allow_pickle=True)
        elif files.endswith(filename + '_opto_response_' + sex + '_side_vel' + '.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response_side_velocity = np.load(path, allow_pickle=True)
        elif files.endswith(filename + '_opto_response_' + sex + '_fwd_vel' + '.npy'):
            path = r'{}/{}'.format(Dir, files)
            opto_response_fwd_velocity = np.load(path, allow_pickle=True)
    if 'monocular' in filename:
        print(opto_response.shape)
        opto_response, opto_response_speed, opto_response_angular_speed, opto_response_side_velocity, opto_response_fwd_velocity = \
            remove_bad_trials([opto_response, opto_response_speed, opto_response_angular_speed, opto_response_side_velocity, opto_response_fwd_velocity], p=p, strict=True)
    else:
        opto_response, opto_response_speed, opto_response_angular_speed, opto_response_side_velocity, opto_response_fwd_velocity = \
            remove_bad_trials([opto_response, opto_response_speed, opto_response_angular_speed, opto_response_side_velocity, opto_response_fwd_velocity], p=p)
    time_periods = 3
    return optomotor_basic_analysis(opto_response, opto_response_speed, opto_response_angular_speed, opto_response_side_velocity, opto_response_fwd_velocity, Dir, filename, sex, time_periods, p=p,
                             condn_indices=condn_indices, ignore_values=ignore_values, variables=variables)


def optomotor_basic_analysis(opto_response, opto_response_speed, opto_response_angular_speed, opto_response_side_velocity, opto_response_fwd_velocity, Dir, filename, p=2, condn_indices=[1, 0], ignore_values = None, variables=['stim_left', 'stim_right']):
    Dir_processed_data = os.path.join(Dir, 'ProcessedData')
    if not os.path.isdir(Dir_processed_data):
        os.mkdir(Dir_processed_data)
    Dir_speed = os.path.join(Dir, 'Speed2')
    if not os.path.isdir(Dir_speed):
        os.mkdir(Dir_speed)
    Dir_OMI = os.path.join(Dir, 'OMI')
    if not os.path.isdir(Dir_OMI):
        os.mkdir(Dir_OMI)
    Dir_basic = os.path.join(Dir, 'Basic')
    if not os.path.isdir(Dir_basic):
        os.mkdir(Dir_basic)
    ##
    indices, rejected_indices = reject_nonresponsive_trials(opto_response_angular_speed, opto_response_speed, p=p, stim_period=stim_period)
    print('{bad} trials rejected from a total {total} ({percent}) for analysis'.format(bad=opto_response.shape[0]-indices.shape[0], total=opto_response.shape[0],
                                                                                    percent=(opto_response.shape[0]-indices.shape[0])/opto_response.shape[0]))
    opto_response = opto_response[indices]
    opto_response_speed = opto_response_speed[indices]
    opto_response_angular_speed = opto_response_angular_speed[indices]
    opto_response_side_velocity = opto_response_side_velocity[indices]
    opto_response_fwd_velocity = opto_response_fwd_velocity[indices]

    data, conditions, indices = split_data_by_conditions([opto_response], [opto_response[:, x] for x in condn_indices], ignore_values=ignore_values)

    indices, conditions = indices_by_condition([opto_response[:, x] for x in condn_indices], ignore_values=ignore_values)
    opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=opto_response[:, p+4:])
    ##
    # fly_id = opto_response[:, p + 2]
    # fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
    # opto_response_list_new = [opto_response_list[i] for i in indices]
    # np.save(os.path.join(Dir, 'data.npy'), opto_response_list)
    # opto_response_angular_speed_list_new = [opto_response_angular_speed[i] for i in indices]
    #
    # draw_trace_each_fly(opto_response_angular_speed_list_new, fly_id_indices, Dir_basic, filename + 'angular_speed', make_labels(conditions),
    #                     top=150, bottom=-150, cmap='nipy_spectral', smoothed='Y', ls=[], sem_or_not='N')
    # draw_trace_each_fly(opto_response_list_new, fly_id_indices, Dir_basic, filename + 'cumulative_turns', make_labels(conditions),
    #                     top=500, bottom=-500, cmap='nipy_spectral', smoothed='Y', ls=[], sem_or_not='N')

    # ## plot for different stimulus conditions without including information of expected direction
    # fig, ax, opto_response_mean, opto_response_sem, opto_response_median, angular_speed_mean, angular_speed_sem, angular_speed_median = draw_optomotor_response(opto_response_list,
    #                 np.array(opto_response_angular_speed), np.clip(opto_response_speed[:, p+2:], 0, 100), Dir_basic, indices, filename=filename, labels=arrow_labels(conditions), top=1000,
    #                 bottom=-1000, cmap=newcmp, ls= ls, rotate='Y', stim_period=stim_period)

    ## add the expected response collumn after p collumns, defaut is 1
    opto_response, exp_type, labels = get_expected_response_hemi(opto_response, p=p, translation='N')

    ## multiply by expected direction such that expected values are positive and opposite responses are negative
    opto_response_list = np.multiply(opto_response_list, opto_response[:, p].reshape(-1, 1) * -1)
    opto_response_angular_speed_list = np.multiply(opto_response_angular_speed[:, p+3:], opto_response[:, p].reshape(-1, 1) * -1)
    opto_response_side_velocity_list = np.multiply(opto_response_side_velocity[:, p+3:], opto_response[:, p].reshape(-1, 1) * -1)
    opto_response_fwd_velocity_list = opto_response_fwd_velocity[:, p+3:]

    ## cluster by BtF, FtB and full rotation and plot separtely
    indices, conditions = indices_by_condition([exp_type], ignore_values=[[5]]) ## take FtB, BtF and full_rotation experiments
    n_conditions = conditions[0][0].shape[0]
    if n_conditions != 0:
        # draw_trace_plot([opto_response_list[i] for i in indices], Dir, filename + 'final_opto_response', labels, top=500, bottom=-500, draw_trace='N', cmap=ListedColormap(color*n_conditions), smoothed='N',
        #                 highlight=None, ls=['solid', 'solid', 'solid'], sem_or_not='Y', normalize='N', rotate='N', sort=[0], draw_median='N', start=180, end=-1, axvline=120,
        #                 PDF='SVG', filter='N')
        # draw_trace_plot([opto_response_angular_speed_list[index] for index in indices], Dir_basic, filename + 'angular_speed_variance', labels=labels, top=300, bottom=0,
        #                 cmap=ListedColormap(color*n_conditions), smoothed='Y', highlight=[0, 0], ls=['dotted', 'dashed', 'solid'], draw_trace='N', draw_linear_sum='N',
        #                  draw_median = 'N', draw_mean='N', draw_var='Y', start=180, end=-1, axvline=120, PDF='N')
        # fig, ax, opto_response_mean_hemi, opto_response_sem_hemi, opto_response_median_hemi, angular_speed_mean_hemi, angular_speed_sem_hemi, angular_speed_median_hemi = \
        #     draw_optomotor_response(opto_response_list, opto_response_angular_speed_list, np.clip(opto_response_speed[:, p+2:], 0, 100), Dir_basic, indices,
        #                             filename=filename+'dir_corrected', labels=labels, top=500, bottom=-500, cmap=ListedColormap(color*n_conditions),
        #                             ls=['dotted', 'dashed', 'solid', 'solid', 'solid', 'solid'], rotate='N', stim_period=stim_period, PDF='Y')
        data = [[opto_response_speed[i][:, p+3:] for i in indices]]
        total_data = {'variables': variables, 'conditions': conditions, 'angular_speed': np.array([angular_speed_mean_hemi, angular_speed_sem_hemi, angular_speed_median_hemi]),
                      'opto_response': np.array([opto_response_mean_hemi, opto_response_sem_hemi, opto_response_median_hemi]),
                      'angular_speed_list': [np.mean(opto_response_angular_speed_list[index][:, stim_period:stim_period + stim_period], axis=1) for index in indices],
                      'opto_response_list': [opto_response_list[index][:, stim_period + stim_period] for index in indices]}
        np.save(os.path.join(Dir_processed_data, filename + '_total_data_hemi.npy'), total_data, allow_pickle=True)

        fly_id = opto_response[:, p + 2]
        fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
        opto_response_list_new = [opto_response_list[i] for i in indices]
        opto_response_angular_speed_list_new = [opto_response_angular_speed_list[i] for i in indices]
        opto_response_side_velocity_list_new = [opto_response_side_velocity_list[i] for i in indices]
        opto_response_fwd_velocity_list_new = [opto_response_fwd_velocity_list[i] for i in indices]
        _, _, _, _, _, cumulative_error_list, prediction_error = draw_trace_plot(opto_response_angular_speed_list_new, Dir_basic, filename + 'angular_speed_linear_comparison',
                        labels=labels, top=150, bottom=-50, cmap=ListedColormap(color*n_conditions), smoothed='Y', highlight=[0, 0], ls=['dotted', 'dashed', 'solid']*2, draw_trace='N',
                        draw_linear_sum='Y', draw_median = 'N', start=180, end=-1, axvline=120, fly_indices=fly_id_indices, PDF='Y')
        draw_trace_polar_plot([opto_response_list[i] for i in indices], Dir, filename, labels, fly_indices=fly_id_indices)

        draw_rastor_plot(opto_response_angular_speed_list_new[::-1], Dir_basic, filename + '_angular_speed_list_unsorted_eachflyPiYG5', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
                         vmin=-150, vmax=150, labels=labels[::-1], smoothed='Y', normalize='N', sort=[], trace_average='Y', start=180, end=-1, axvline=120, fly_indices=fly_id_indices[::-1], PDF='SVG')
        draw_rastor_plot(opto_response_side_velocity_list_new[::-1], Dir_basic, filename + '_side_velocity_eachflyPiYG5', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
                         vmin=-10, vmax=10, labels=labels[::-1], smoothed='Y', normalize='N', sort=[], trace_average='Y', start=180, end=-1, axvline=120, fly_indices=fly_id_indices[::-1], PDF='SVG')
        draw_rastor_plot(opto_response_fwd_velocity_list_new[::-1], Dir_basic, filename + '_fwd_velocity_eachflyPiYG5', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
                         vmin=-50, vmax=50, labels=labels[::-1], smoothed='Y', normalize='N', sort=[], trace_average='Y', start=180, end=-1, axvline=120, fly_indices=fly_id_indices[::-1], PDF='SVG')

        ##new addition
        # Dir = os.path.join(r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures', 'multiple_clusters', filename)
        # if not os.path.isdir(Dir):
        #     os.mkdir(Dir)
        # split_and_plot_rastor(opto_response_angular_speed_list_new[::-1], opto_response_list_new[::-1],  fly_id_indices[::-1], Dir, filename+"_manual_clustered", labels[::-1], 'SVG')
        # return opto_response_angular_speed_list_new, fly_id_indices

        fly_id = opto_response[:, p + 2]
        fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
        opto_response_angular_speed_list_new = [opto_response_angular_speed_list[i] for i in indices]
        draw_rastor_plot(opto_response_angular_speed_list_new[::-1], Dir_basic, filename + '_angular_speed_list_random', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
                         vmin=-150, vmax=150, labels=labels[::-1], smoothed='Y', normalize='N', sort=[0], trace_average='Y', start=180, end=-1, axvline=120, fly_indices=fly_id_indices[::-1], PDF='SVG', to_plot=-200)
        fly_id = opto_response[:, p + 2]
        fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
        opto_response_angular_speed_list_new = [opto_response_angular_speed_list[i] for i in indices]
        draw_rastor_plot(opto_response_angular_speed_list_new[::-1], Dir_basic, filename + '_angular_speed_list_random2', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
                         vmin=-150, vmax=150, labels=labels[::-1], smoothed='Y', normalize='N', sort=[0], trace_average='Y', start=180, end=-1, axvline=120, fly_indices=fly_id_indices[::-1], PDF='SVG', to_plot=-200)

        total_data = {'cumulative_error_list': cumulative_error_list, 'prediction_error': prediction_error}
        np.save(os.path.join(Dir_processed_data, filename + '_cumulative_error_list.npy'), total_data, allow_pickle=True)

        # doing this again because fly_indices is changed by the previous function
        fly_id = opto_response[:, p + 2]
        fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
        opto_response_list_new = [opto_response_list[i] for i in indices]
        opto_response_angular_speed_list_new = [opto_response_angular_speed_list[i] for i in indices]
        #
        draw_trace_each_fly(opto_response_list_new, fly_id_indices, Dir_basic, filename, labels, top=750, bottom=-750, cmap='tab20', smoothed='N', ls=[], turn_or_angular_speed=True,
                            axvline=120, start=180, end=-1)

        draw_rastor_plot_with_cluster_id(opto_response_angular_speed_list_new, Dir_basic, filename+'angular_speed_rastor_dly_ids', cmap=ListedColormap(all_color_data['PiYG5']), vmin=-200, vmax=200, labels=labels,
                     smoothed='Y', normalize= 'N',  PDF='SVG', axvline=300, cluster_ids=[fly_id_indices], cluster_cmap='tab20b')
    return 1
