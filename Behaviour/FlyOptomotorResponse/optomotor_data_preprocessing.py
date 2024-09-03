import os
import numpy as np
from scipy.signal import savgol_filter
from matplotlib.colors import ListedColormap

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Import your custom modules here
import read_folder
from detect_saccades import detect_saccades_cwt
from check_saccade import check_saccade
from stim_triggered_average import stim_triggered_average
from frame_timeperiod import frame_timeperiod
from stimulus_pattern_hemi import stimulus_pattern_hemi
from indices_by_condition import indices_by_condition
from draw_position_heatmap import draw_position_heatmap
from get_expected_response_hemi import get_expected_response_hemi


def main(Dir, sex='Y', PARAMETERS = ['frequency', 'bars','stimulus'] , fly_id=None, saccade='Y', normal='Y', new=True):
    plt.interactive(False)
    flystack = []
    Dir = Dir
    flystack = read_folder.search_data(Dir, flystack, find='corrected', ignore='speed', fix_csv='N')
    name = os.path.split(Dir)[-1]
    print(name)
    if fly_id == None:
        fly_id = name.split('_')[-1]
    else:
        pass

    for j in range(len(flystack)):
        flydata = flystack[j]
        fly_sex = flydata.fly['sex']
        if sex == fly_sex:
            pass
        elif sex == 'Y':
            pass
        else:
            return [0]
        start = 10
        for m in range(len(flydata.track)):
            pos_x = np.array(flydata.track[m]['pos_x'])
            pos_y = np.array(flydata.track[m]['pos_y'])

            fig, ax = plt.subplots()
            ax.plot(pos_x[start * 60 * 60:int(2+start) * 60 * 60], pos_y[start * 60 * 60:int(2+start) * 60 * 60], lw=1, color='k')
            circle2 = mpatches.Circle((500, 500), radius=475, color='gray', fill=False, lw=2)
            ax.add_artist(circle2)
            ax.set_frame_on(False)
            ax.set_xlim(left=0, right=1000)
            ax.set_ylim(bottom=0, top=1000)
            ax.set(aspect=1)
            ax.set_xticks([])
            ax.set_yticks([])
            fig.savefig(os.path.join(Dir, 'trajectory.svg'), dpi=600, transparent=True)
            fig.savefig(os.path.join(Dir, 'trajectory.png'))

            draw_position_heatmap(pos_x, pos_y, Dir=Dir, bins_rad=50, bins_theta=50, radius=512, file='position_heatmap_male.png', color_map='afmhot', vmin=0, vmax=0.00001)
            ori = np.array(flydata.track[m]['ori'])
            ori_360 = np.remainder(ori, 360)

            direction_left = np.array(flydata.track[m]['direction_left'])
            direction_right = np.array(flydata.track[m]['direction_right'])
            try:
                direction_center = np.array(flydata.track[m]['direction_center'])
                print('the stimulus has 3 sections')
                sections = 3
            except:
                direction_center = np.zeros(direction_left.shape[0])
                print('the stimulus has 2 sections')
                sections = 2

            try:
                direction_rear = np.array(flydata.track[m]['direction_rear'])
                print('the stimulus has 4 sections')
                sections = 4
            except:
                direction_rear = np.zeros(direction_left.shape[0])
                print('the stimulus has 3 sections')
                sections = 3

            trial_period = stim_period * 4
            no_of_trials = len(direction_left) // trial_period
            print('number of trirlas {trials}'.format(trials=no_of_trials))
            # frame_changes = np.arange(stim_period, len(direction_left), trial_period//2)
            ##
            frame_changes = stimulus_pattern_hemi(direction_left, direction_right, direction_center, direction_rear)
            frame_changes = np.array(frame_changes)
            ##
            ## take values after 60 frames
            stim_change_left = np.array([np.median(direction_left[x:x+60]) for x in frame_changes]).astype(np.int)
            stim_change_right = np.array([np.median(direction_right[x:x+60]) for x in frame_changes]).astype(np.int)
            if sections == 3:
                stim_change_center = np.array([np.median(direction_center[x:x + 60]) for x in frame_changes]).astype(np.int)
            elif sections == 4:
                stim_change_center = np.array([np.median(direction_center[x:x + 60]) for x in frame_changes]).astype(np.int)
                stim_change_rear = np.array([np.median(direction_rear[x:x + 60]) for x in frame_changes]).astype(np.int)
            else:
                pass

            try:
                fly_lost = np.array(flydata.track[m]['fly_lost'])
                correct_ori = np.array(flydata.track[m]['correct_ori'])
                ## new method to fix bad trials
                # fix the frames with wrong ori by switching left and right
                # but, only in trials where the direction was flipped throughout
                ## find trials that were flipped throughout
                for i, frames in enumerate(frame_changes):
                    if np.argwhere(fly_lost[frames + 1:frames + stim_period] == -1).flatten().shape[0] >= 285:
                    # if np.all(fly_lost[frames+1:frames + stim_period] == -1):
                        correct_right_ori = copy.deepcopy(direction_left[frames:frames + stim_period])
                        direction_left[frames:frames + stim_period] = direction_right[frames:frames + stim_period]
                        direction_right[frames:frames + stim_period] = correct_right_ori
                        #
                        if sections==4 or sections==3:
                            correct_center_ori = copy.deepcopy(direction_rear[frames:frames + stim_period])
                            direction_rear[frames:frames + stim_period] = direction_center[frames:frames + stim_period]
                            direction_center[frames:frames + stim_period] = correct_center_ori
                        #
                        if 'new' in name:
                            fly_lost[frames:frames + stim_period] = 0  ## fix the fly_lost data, not incorrect anymore
                        else:
                            fly_lost[frames:frames + stim_period] = 5  ## fix the fly_lost data, 5 means that the data was wrong but is fixed now
                # wrong_ori_frames = np.argwhere(fly_lost == -1)  ## the stimulus was presented opposite to what it should be
                # fly_lost_frames = np.argwhere(fly_lost == 1)  ## the fly was lost due to jumping and hitting the wall
                # ## derive the trials where orientation is determined incorrectly
                # incorrect_trials = wrong_ori_frames // stim_period
                # ## array of trials with incorrect ori_frames
                # incorrect_trials_set = np.unique(incorrect_trials)
                ## bad_frames where the fly is lost or the orientation was incorrectly determined
                bad_frames = np.argwhere((fly_lost == 1) | (fly_lost == -1))
                fixed_frames = np.argwhere(fly_lost == 5)
                ## derive the trials where there were bad frames
                bad_trials = bad_frames // stim_period
                fixed_trials = fixed_frames // stim_period
                ## array of trials with bad frames
                bad_trials_set, bad_trials_count = np.unique(bad_trials, return_counts=True)
                bad_trials_set = bad_trials_set[np.argwhere(bad_trials_count>15).flatten()]
                fixed_trials_set, fixed_trials_count = np.unique(fixed_trials, return_counts=True)
                fixed_trials_set = fixed_trials_set[np.argwhere(fixed_trials_count>15).flatten()]
                
                ## assign good and bad periods to stimulus and pre-stimulus periods
                bad_trial_stim = []
                bad_trial_non_stim = []
                fixed_trial_stim = []
                fixed_trial_non_stim = []
                for i in bad_trials_set:
                    count = np.count_nonzero(bad_trials == i)
                    if count > 10:
                        if i % 2 != 0: ## even half_trials are periods with rotation
                            bad_trial_stim.append(i // 2)
                        else: ## odd half_trials are periods without rotation
                            bad_trial_non_stim.append(i // 2)
                for i in fixed_trials_set:
                    count = np.count_nonzero(fixed_trials == i)
                    if count > 10:
                        if i % 2 != 0: ## even half_trials are periods with rotation
                            fixed_trial_stim.append(i // 2)
                        else: ## odd half_trials are periods without rotation
                            fixed_trial_non_stim.append(i // 2)
            except:
                pass

            ## take values after 60 frames
            stim_change_left = np.array([np.median(direction_left[x:x+60]) for x in frame_changes]).astype(np.int)
            stim_change_right = np.array([np.median(direction_right[x:x+60]) for x in frame_changes]).astype(np.int)
            if sections == 3:
                stim_change_center = np.array([np.median(direction_center[x:x + 60]) for x in frame_changes]).astype(np.int)
            elif sections == 4:
                stim_change_center = np.array([np.median(direction_center[x:x + 60]) for x in frame_changes]).astype(np.int)
                stim_change_rear = np.array([np.median(direction_rear[x:x + 60]) for x in frame_changes]).astype(np.int)
            else:
                pass

            timestamp = np.array(flydata.track[m]['timestamp'])
            timeperiod = frame_timeperiod(timestamp)
            # return timeperiod
            try:
                try:
                    freq_left = stim_params_list_json(flydata, 'temp', 0)[0] ## speed
                    bars_left = stim_params_list_json(flydata, 'spatial', 0)[0] ## number of dots
                    contrast_left = stim_params_list_json(flydata, 'contrast', 0)[0] ## contrast
                    optic_flow_left = (360/np.array(bars_left))*np.array(freq_left)
                    total_optic_flow_left = optic_flow_left * 10

                    freq_right = stim_params_list_json(flydata, 'temp', 0)[1]
                    bars_right = stim_params_list_json(flydata, 'spatial', 0)[1]
                    contrast_right = stim_params_list_json(flydata, 'contrast', 0)[1]
                    optic_flow_right = (360/np.array(bars_right))*np.array(freq_right)
                    total_optic_flow_right = optic_flow_right * 10
                except:
                    freq_left = list(np.array(stim_params_list_json(flydata, 'temporal_frequency', 0))[:, 0]) ## speed
                    bars_left =  list(np.array(stim_params_list_json(flydata, 'spatial_frequency', 0))[:, 0]) ## number of dots
                    contrast_left =  list(np.array(stim_params_list_json(flydata, 'contrast', 0))[:, 0]) ## contrast
                    optic_flow_left = (360/np.array(bars_left))*np.array(freq_left)
                    total_optic_flow_left = optic_flow_left * 10

                    freq_right = list(np.array(stim_params_list_json(flydata, 'temporal_frequency', 0))[:, 1]) ## speed
                    bars_right =  list(np.array(stim_params_list_json(flydata, 'spatial_frequency', 0))[:, 1]) ## number of dots
                    contrast_right =  list(np.array(stim_params_list_json(flydata, 'contrast', 0))[:, 1]) ## contrast
                    optic_flow_right = (360/np.array(bars_right))*np.array(freq_right)
                    total_optic_flow_right = optic_flow_right * 10
            except:
                rod = stim_params_list_json(flydata, 'rod', 0)  ## speed
                N_dots = stim_params_list_json(flydata, 'N', 0)  ## number of dots
                freq_left = stim_params_list_json(flydata, 'angle', 0)
                contrast_right = [0] * len(freq_left)
                contrast_left = [0] * len(freq_left)
                freq_right = [0] * len(freq_left)
                total_optic_flow_left = [0] * len(freq_left)
                total_optic_flow_right = [0] * len(freq_left)

            try:
                closedloop_gain = np.array(stim_params_list_json(flydata, 'stim_attributes', 0)[0])
                closedloop_gain = closedloop_gain[:, 1]
                closedloop_gain = closedloop_gain.astype(np.float)
                closedloop_gain = np.multiply(closedloop_gain, 100)
            except:
                print('There is no closedloop gain')
            try:
                stimulus_bar = np.array(stim_params_list_json(flydata, 'stim', 0)[0], dtype='')[:, 3]
            except:
                print('nothing found')
                stimulus_bar = np.zeros(len(freq_left)).astype(int)

            # stim_type = np.array(stim_params_list_json(flydata, 'stim', 0)[0])
            # stim_type = stim_type[:, 0]
            # if stim_type[0] == 'RDK':
            #     ## for random dots, freq==speed and bars==density of dots
            #     pass
            # elif stim_type[0] == 'ring_pinwheel':
            #     ring_size = np.array(stim_params_list_json(flydata, 'stim', 0)[0])
            #     ring_size = ring_size[:,2]
            # else:
            #     pass
            # ##########################
            # stimulus_type = []
            # for i in range(stim_type.shape[0]):
            #     if stim_type[i] == 'sine_pinwheel':
            #         stimulus_type.append(0)
            #     elif 'sine_pinwheel_new' in stim_type[i]:
            #         stimulus_type.append(1)
            stimulus_type = [1] * len(freq_left)
            ##########################
            try:
                color = np.array(stim_params_list_json(flydata, 'stim', m)[0])[:, 2]
                color_set = list(np.unique(color))
                color_set.sort()
                color_stim = []
                for i in color:
                    for j in range(len(color_set)):
                        if i == color_set[j]:
                            color_stim.append(j)
                            break
            except:
                print('Background color not given')
            ##########################
            #############################
            # for i in range(len(rod)):
            #     if rod[i] == 20:
            #         # print('red')
            #         rod[i] = 25
            #     if rod[i] == 2:
            #         # print('red')
            #         rod[i] = 1
            # ############################
            ori_smooth, ori_360_smooth, angular_speed, speed = compute_velocities(timeperiod, pos_x, pos_y, ori, no_of_trials, trial_period)
            heading_direction, dist = locomotion_direction(pos_x, pos_y, correct_ori, 5)
            fwd_velocity, side_velocity = speed_components(speed, heading_direction, correct_ori)

            # return angular_speed, speed, heading_direction, correct_ori
            # check_tracking_results(Dir, correct_ori, heading_direction) #comment out for analysis

            acceleration = np.diff(speed[1:])
            acceleration = velocity(acceleration, timeperiod[2:])
            acceleration = np.insert(acceleration, 0, 0)
            acceleration = np.insert(acceleration, 0, 0)
            #########################
            r, theta = cart2pol(np.subtract(pos_x, 512), np.subtract(pos_y, 512))
            filtered_angular_speed = copy.deepcopy(angular_speed)
            i = 0
            while i < len(filtered_angular_speed) - 30:
                if r[i] > 300:
                    filtered_angular_speed[i:i + 30] = np.zeros(30)
                    i += 30
                else:
                    pass
                i += 1
            filtered_ori_smooth = copy.deepcopy(ori_smooth)
            i = 0
            while i < len(filtered_angular_speed) - 30:
                if r[i] > 300:
                    filtered_ori_smooth[i:i + 30] = np.zeros(30)
                    i += 30
                else:
                    pass
                i += 1

            # lost_fly_frames = np.where((pos_x==0) & (pos_y==0))[0]
            # if lost_fly_frames.shape[0]==0:
            #     n = len(freq_left)
            # else:
            #     n = np.amin(lost_fly_frames)//1200
            if new:
                n = len(freq_left)
                print(n, 'n')
            else:
                n = 2*len(freq_left)
            stim_change_left = stim_change_left[:n]
            stim_change_right = stim_change_right[: n]
            if sections == 3:
                stim_change_center = stim_change_center[:n]
                sections = 4
                stim_change_rear = np.zeros((1, n))
            else:
                pass
            if sections == 4:
                stim_change_center = stim_change_center[:n]
                stim_change_rear = stim_change_rear[:n]
            else:
                pass
            stim_change = np.ones(stim_change_right.shape[0])
            ##
            trial_quality = np.zeros(n)
            bad_trial_non_stim = np.array(bad_trial_non_stim)
            bad_trial_stim = np.array(bad_trial_stim)
            bad_trial_non_stim = bad_trial_non_stim[np.array(bad_trial_non_stim) < n].flatten()
            bad_trial_stim = bad_trial_stim[np.array(bad_trial_stim) < n].flatten()
            ##
            if bad_trial_non_stim.shape[0]!=0:
                trial_quality[bad_trial_non_stim] = -1 #means that only the non_stim period is bad
            if bad_trial_stim.shape[0]!=0:
                trial_quality[bad_trial_stim] = 1 #means that the stim period is bad, the non_stim period might be bad too
            ##
            fixed_trial_non_stim = np.array(fixed_trial_non_stim)
            fixed_trial_stim = np.array(fixed_trial_stim)
            fixed_trial_non_stim = fixed_trial_non_stim[np.array(fixed_trial_non_stim) < n].flatten()
            fixed_trial_stim = fixed_trial_stim[np.array(fixed_trial_stim) < n].flatten()
            ##
            if fixed_trial_non_stim.shape[0]!=0:
                trial_quality[fixed_trial_non_stim] = -5 #means that only the non_stim period is bad
            if fixed_trial_stim.shape[0]!=0:
                trial_quality[fixed_trial_stim] = 5 #means that the stim period is bad, the non_stim period might be bad too
            print(np.argwhere(trial_quality==5).shape, 'here')
            ##
            if sections == 2:
                full_rotation_indices = np.argwhere(stim_change_right == stim_change_left).flatten()
                trial_quality[full_rotation_indices] = 0 ## the orientation of the trial does not matter for full rotation stimulus

            ##
            # stim_change = []
            # for i in range(stim_change_right.shape[0]):
            #     if stim_change_right[i] == stim_change_left[i]:
            #         stim_change.append(stim_change_right[i])
            #     else:
            #         if stim_change_left[i] == 0:
            #             stim_change.append(stim_change_right[i])
            #         elif stim_change_right[i] == 0:
            #             stim_change.append(stim_change_left[i])
            #         else:
            #             stim_change.append(1)
            # stim_change = np.array(stim_change)

            # pos_x_new = pos_x[:n * trial_period]
            # pos_y_new = pos_y[:n * trial_period]
            # pos_x_reshaped = pos_x_new.reshape(pos_x_new.shape[0] // trial_period//2, trial_period//2)
            # pos_y_reshaped = pos_y_new.reshape(pos_y_new.shape[0] // trial_period//2, trial_period//2)
            # pos_reshaped = np.stack((pos_x_reshaped, pos_y_reshaped), axis = 2)
            pos_x_reshaped = stim_triggered_average(frame_changes, pos_x, stim_period, stim_period)
            pos_y_reshaped = stim_triggered_average(frame_changes, pos_y, stim_period, stim_period)

            pos_reshaped = np.stack((pos_x_reshaped, pos_y_reshaped), axis=2)

            if saccade == 'Y':
                if sections == 4:
                    parameter = define_parameters4(stim_change_left, stim_change_right, stim_change_center, stim_change_rear, contrast_left, contrast_right, stimulus_type, freq_left, freq_right, PARAMETERS=PARAMETERS)
                    saccade_information = fly_saccade_information(filtered_angular_speed, filtered_ori_smooth, frame_changes, stim_change, trial_period, parameter, trial_quality)
                elif sections == 3:
                    parameter = define_parameters(stim_change_left, stim_change_right, stim_change_center, contrast_left,contrast_right, stimulus_type, freq_left, freq_right, PARAMETERS=PARAMETERS)
                    saccade_information = fly_saccade_information(filtered_angular_speed, filtered_ori_smooth, frame_changes, stim_change, trial_period, parameter, trial_quality)
                else:
                    parameter = define_parameters2(stim_change_left, stim_change_right, contrast_left,contrast_right, stimulus_type, freq_left, freq_right, PARAMETERS=PARAMETERS)
                    saccade_information = fly_saccade_information(filtered_angular_speed, filtered_ori_smooth, frame_changes, stim_change, trial_period, parameter, trial_quality)
                if normal=='N':
                    return saccade_information, pos_x_reshaped, pos_y_reshaped, speed
                else:
                    pass
            else:
                saccade_information = 1
            pos_reshaped = pos_reshaped[:n]
            orientation = stim_triggered_average(frame_changes, correct_ori, stim_period, stim_period)  #### change to 600 here
            orientation = orientation[:n]
            opto_response = stim_triggered_average(frame_changes, ori_smooth, stim_period, stim_period)  #### change to 600 here
            opto_response = opto_response[:n]
            opto_response_speed = stim_triggered_average(frame_changes, speed, stim_period, stim_period)
            opto_response_speed = opto_response_speed[:n]
            opto_response_angular_speed = stim_triggered_average(frame_changes, angular_speed, stim_period, stim_period)
            opto_response_angular_speed = opto_response_angular_speed[:n]
            opto_response_fwd_vel = stim_triggered_average(frame_changes, fwd_velocity, stim_period, stim_period)
            opto_response_side_vel = stim_triggered_average(frame_changes, side_velocity, stim_period, stim_period)
            opto_response_fwd_vel = opto_response_fwd_vel[:n]
            opto_response_side_vel = opto_response_side_vel[:n]
            ##########################
            if sections==4:
                parameter = define_parameters4(stim_change_left, stim_change_right, stim_change_center, stim_change_rear, contrast_left, contrast_right, stimulus_type, freq_left, freq_right, PARAMETERS=PARAMETERS)
            elif sections==3:
                parameter = define_parameters(stim_change_left, stim_change_right, stim_change_center, contrast_left, contrast_right, stimulus_type, freq_left, freq_right, PARAMETERS=PARAMETERS)
            else:
                parameter = define_parameters2(stim_change_left, stim_change_right, contrast_left, contrast_right, stimulus_type, freq_left, freq_right, PARAMETERS=PARAMETERS)
            p=0
            ## insert the parameters
            for i in range(len(parameter)):
                opto_response = np.insert(opto_response, i, parameter[i], 1)
                opto_response_angular_speed = np.insert(opto_response_angular_speed, i, parameter[i], 1)
                opto_response_speed = np.insert(opto_response_speed, i, parameter[i], 1)
                opto_response_fwd_vel = np.insert(opto_response_fwd_vel, i, parameter[i], 1)
                opto_response_side_vel = np.insert(opto_response_side_vel, i, parameter[i], 1)
                orientation = np.insert(orientation, i, parameter[i], 1)
                p=i
            print(p, 'p')
            ## insert the expected direction of response, default is 1
            opto_response = np.insert(opto_response, p+1, np.array(stim_change), 1)
            orientation = np.insert(orientation, p + 1, np.array(stim_change), 1)
            opto_response_angular_speed = np.insert(opto_response_angular_speed, p+1, np.array(stim_change), 1)
            opto_response_speed = np.insert(opto_response_speed, p+1, np.array(stim_change), 1)
            opto_response_fwd_vel = np.insert(opto_response_fwd_vel, p + 1, np.array(stim_change), 1)
            opto_response_side_vel = np.insert(opto_response_side_vel, p + 1, np.array(stim_change), 1)
            opto_response = np.insert(opto_response, p+2, np.repeat(total_optic_flow_left[:n], 1), 1)
            #### Insert the fly_id
            opto_response = np.insert(opto_response, p+3, np.ones(opto_response.shape[0])*np.asarray(fly_id, dtype='float64'), 1)
            orientation = np.insert(orientation, p + 3, np.ones(orientation.shape[0]) * np.asarray(fly_id, dtype='float64'), 1)
            opto_response_speed = np.insert(opto_response_speed, p+2, np.ones(opto_response_speed.shape[0])*np.asarray(fly_id, dtype='float64'), 1)
            opto_response_angular_speed = np.insert(opto_response_angular_speed, p+2, np.ones(opto_response_angular_speed.shape[0])*np.asarray(fly_id, dtype='float64'), 1)
            opto_response_side_vel = np.insert(opto_response_side_vel, p + 2, np.ones(opto_response_side_vel.shape[0]) * np.asarray(fly_id, dtype='float64'), 1)
            opto_response_fwd_vel = np.insert(opto_response_fwd_vel, p + 2, np.ones(opto_response_fwd_vel.shape[0]) * np.asarray(fly_id, dtype='float64'), 1)
            ### insert trial quality 0=perfect, 1=stim_phase bad, -1=non_stim_phase bad, 5=stim phase bad but fixed, -5=non_stim_phase bad but fixed
            print(np.argwhere(trial_quality == 5).shape, 'here')
            opto_response = np.insert(opto_response, p+4, trial_quality, 1)
            orientation = np.insert(orientation, p+4, trial_quality, 1)
            opto_response_speed = np.insert(opto_response_speed, p+3, trial_quality, 1)
            opto_response_angular_speed = np.insert(opto_response_angular_speed, p+3, trial_quality, 1)
            opto_response_side_vel = np.insert(opto_response_side_vel, p+3, trial_quality, 1)
            opto_response_fwd_vel = np.insert(opto_response_fwd_vel, p+3, trial_quality, 1)
            # pos_reshaped = np.insert(pos_reshaped, 0, trial_quality, 1)
            ###
            opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0]//2]), axis=1, arr=opto_response[:, 5:])
            indices, conditions = indices_by_condition([opto_response[:, 0]]) ## change p to 0
        print('Done.........')
        print(p)
        opto_response, exp_type, labels = get_expected_response_hemi(opto_response, p=p+1)
        print(labels)
        ## multiply by expected direction such that expected values are positive and opposite responses are negative
        opto_response_list = np.multiply(opto_response_list, opto_response[:, p].reshape(-1, 1) * -1)
        opto_response_angular_speed_list = np.multiply(opto_response_angular_speed[:, p + 3:], opto_response[:, p].reshape(-1, 1) * -1)
        ## plot cw and ccw full rotation response
        # fig, ax = plt.subplots(figsize=(3, 3))
        # cw_indices = np.argwhere((opto_response[:, 0]==1) & (opto_response[:, 1]==1)).flatten()
        # ccw_indices = np.argwhere((opto_response[:, 0]==-1) & (opto_response[:, 1]==-1)).flatten()
        # cmap = ListedColormap(all_color_data['RdBu3'])
        # colors = [cmap(0.9), cmap(0.1)]
        # # ax.plot(savgol_filter(np.mean(opto_response_angular_speed_list[cw_indices], axis=0)*-1, 31, 2), c=colors[0], lw=2.5)
        # # ax.plot(savgol_filter(np.mean(opto_response_angular_speed_list[ccw_indices], axis=0), 31, 2), c=colors[1], lw=2.5)
        # ax.set_ylim(top=400, bottom=-400)
        # ax.axvline(x=300, lw=1, ls='--', c='grey')
        # ax.set_yticks([-400, 0, 400])
        # ax.spines['bottom'].set_bounds(high=600, low=0)
        # ax.spines['left'].set_bounds(high=400, low=-400)
        # # fig.savefig(os.path.join(Dir, 'cw_ccw_comparison.svg'), dpi=600, transparent=True)
        # fig.savefig(os.path.join(Dir, 'cw_ccw_comparison.png'))
        # ## cluster by BtF, FtB and full rotation and plot separtely
        # # color = [color_dict['CantonS'][1]]
        # indices, conditions = indices_by_condition([exp_type], ignore_values=[[5]])  ## take FtB, BtF and full_rotation experiments
        # n_conditions = conditions[0][0].shape[0]
        # # if n_conditions != 0:
        # #     color = ['k']
        # #     draw_optomotor_response(opto_response_list, opto_response_angular_speed_list, np.clip(opto_response_speed[:, p + 2:], 0, 100), Dir, indices,
        # #                                 filename=name, labels=labels, top=500, bottom=-500, cmap=ListedColormap(color * n_conditions),
        # #                                 ls=['dotted', 'dashed', 'solid', 'solid', 'solid', 'solid'], rotate='N', stim_period=stim_period, PDF='N')
        # save filename and assigned fly_id
        newdict = {'filename':Dir, 'flyid':fly_id}
        return opto_response, opto_response_speed, opto_response_angular_speed, saccade_information, opto_response_fwd_vel, opto_response_side_vel, pos_reshaped, orientation, newdict
    return flystack


def fly_saccade_information(angular_speed, ori_smooth, frame_changes, stim_change, trial_period, var, trial_quality):
    print(stim_change.shape)
    peaks, peak_values, width_heights, peak_widths = detect_saccades_cwt(angular_speed, 200)
    ## peaks = indices of peak in angular speed
    ## peak_values = angular speed at the peak
    saccade_total_turn, saccade_width_turn = check_saccade(peaks, peak_widths, ori_smooth)
    saccade_turn_all = np.zeros(angular_speed.shape[0])
    saccade_peak_all = np.zeros(angular_speed.shape[0])
    saccade_peak_all[peaks] = peak_values
    saccade_turn_all[peaks] = saccade_total_turn
    n = len(var[0])
    stim_change = stim_change[:n]
    saccade_peak_list = stim_triggered_average(frame_changes, saccade_peak_all, stim_period, stim_period)
    saccade_turn_list = stim_triggered_average(frame_changes, saccade_turn_all, stim_period, stim_period)
    saccade_peak_list = saccade_peak_list[:n, :]
    saccade_turn_list = saccade_turn_list[:n, :]
    saccade_peak_list = np.insert(saccade_peak_list, 0, np.array(stim_change), 1)
    saccade_turn_list = np.insert(saccade_turn_list, 0, np.array(stim_change), 1)
    for i in range(len(var)):
        saccade_peak_list = np.insert(saccade_peak_list, i + 1, var[i], 1)
        saccade_turn_list = np.insert(saccade_turn_list, i + 1, var[i], 1)
    saccade_peak_list = np.insert(saccade_peak_list, i+2, trial_quality, 1)
    saccade_turn_list = np.insert(saccade_turn_list, i+2, trial_quality, 1)
    return saccade_peak_list, saccade_turn_list