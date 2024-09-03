import numpy as np
import math
import copy
import pandas as pd

def velocity(diff_data, time_periods, factor = 1e6):
    """
    calculate velocity from the difference in data and time period
    """
    if (type(diff_data) == list) & (type(time_periods) == list):
        assert(len(diff_data) == len(time_periods)),"Length of the two lists not the same"
    else:
        print(type(diff_data), type(time_periods))
        assert(diff_data.shape[0] == len(time_periods)),"Length of the two lists not the same {}, {}".format(diff_data.shape[0], len(time_periods))

    return np.multiply(np.divide(diff_data, time_periods),factor)


def compute_velocities_courtship_150(timeperiod, pos_x_m, pos_y_m, pos_x_f, pos_y_f, ori_360_m, ori_360_f, 
                                 timefactor=1e6, window=5):
    turns_m = np.subtract((ori_360_m[2:]), (ori_360_m[:-2]))
    turns_m = np.pad(turns_m, [1,1], 'constant')
    for i in range(turns_m.shape[0]):  # to fix angles where the fly completes 360 turns
        if turns_m[i] > 270:
            turns_m[i] = turns_m[i] - 360
        elif turns_m[i] < -270:
            turns_m[i] = 360 + turns_m[i]

    errors = [i for i in range(turns_m.shape[0]) if turns_m[i] > 90 or turns_m[i] < -90]  # frames where the male is making larger than 90 degree turns
    error_diff = [errors[i + 1] - errors[i] for i in range(0, len(errors) - 1, 2)]  # check every alternative incorrect frame to determine the duration of the incorrect bout
    # the logic here is that 2 consecutive erroneuous entries indicate a 360 turn back to correct orientation
    min_error = np.sum(error_diff)  # total length of incorrect bouts
    i = 0
    errors_copy = copy.deepcopy(errors)
    min_point = -1
    while i < len(errors):
        errors.remove(errors[i])  # remove one of the incorrect frames
        error_diff = [errors[i + 1] - errors[i] for i in range(0, len(errors) - 1, 2)]
        total_error = np.sum(error_diff)  # recompute error
        if total_error < min_error:  # check if that reduces the total length of incorrect bouts
            # the logic is that very few frames are incorrect, not all, so we are trying to get that
            min_error = total_error
            min_point = i
        errors = copy.deepcopy(errors_copy)
        i = i + 1
    if min_point == -1:
        pass
    else:
        print('min point{}'.format(min_point))
        min_point = errors_copy[min_point]
        #####################
        for i in range(turns_m.shape[0]):
            if i == min_point:
                continue
            else:
                if turns_m[i] > 90:
                    turns_m[i] = turns_m[i] - 180
                elif turns_m[i] < -90:
                    turns_m[i] = 180 + turns_m[i]

    turns_smooth_m = np.median(np.lib.stride_tricks.sliding_window_view(turns_m, (window,)), axis=1)
    if window % 2 == 0:
        turns_smooth_m = np.pad(turns_smooth_m, [window//2 - 1, window//2], 'constant')
    else:
        turns_smooth_m = np.pad(turns_smooth_m, [window//2, window//2], 'constant')
    turns_mean_smooth_m = np.convolve(turns_smooth_m, np.ones(5) / 5, 'same')
    ori_smooth_m = np.cumsum(turns_smooth_m)   
    ori_smooth_m_shuffled = copy.deepcopy(ori_smooth_m)
    np.random.shuffle(ori_smooth_m_shuffled)

    ang_vel_m = velocity(np.array(turns_m), timeperiod, timefactor)
    ang_vel_smooth_m = np.median(np.lib.stride_tricks.sliding_window_view(ang_vel_m, (window,)), axis=1)
    if window % 2 == 0:
        ang_vel_smooth_m = np.pad(ang_vel_smooth_m, [window//2 - 1, window//2], 'constant')
    else:
        ang_vel_smooth_m = np.pad(ang_vel_smooth_m, [window//2, window//2], 'constant')

    #for female
    turns_f = np.subtract((ori_360_f[2:]), (ori_360_f[:-2]))
    turns_f = np.pad(turns_f, [1,1], 'constant')
    for i in range(turns_f.shape[0]):
        if turns_f[i] > 270:
            turns_f[i] = 360 - turns_f[i]
        elif turns_f[i] < -270:
            turns_f[i] = turns_f[i] + 360

    for i in range(turns_f.shape[0]):
        if turns_f[i] > 90:
            turns_f[i] = turns_f[i] - 180
        elif turns_f[i] < -90:
            turns_f[i] = 180 + turns_f[i]
    
    turns_smooth_f = np.median(np.lib.stride_tricks.sliding_window_view(turns_f, (window,)), axis=1)
    if window % 2 == 0:
        turns_smooth_f = np.pad(turns_smooth_f, [window//2 - 1, window//2], 'constant')
    else:
        turns_smooth_f = np.pad(turns_smooth_f, [window//2, window//2], 'constant')
    turns_mean_smooth_f = np.convolve(turns_smooth_f, np.ones(5) / 5, 'same')
    ori_smooth_f = np.cumsum(turns_smooth_f)
    ori_smooth_f_shuffled = copy.deepcopy(ori_smooth_f)
    np.random.shuffle(ori_smooth_f_shuffled)
    
    ang_vel_f = velocity(turns_smooth_f, timeperiod, timefactor)
    ang_vel_smooth_f = np.median(np.lib.stride_tricks.sliding_window_view(ang_vel_f, (window,)), axis=1)
    if window % 2 == 0:
        ang_vel_smooth_f = np.pad(ang_vel_smooth_f, [window//2 - 1, window//2], 'constant')
    else:
        ang_vel_smooth_f = np.pad(ang_vel_smooth_f, [window//2, window//2], 'constant')
    ang_vel_mean_smooth_f = np.convolve(ang_vel_smooth_f, np.ones(5) / 5, 'same')

    desired_female_pos = []
    for i in range(len(pos_x_m)):
        ep1 = (int(pos_x_m[i] + 60 * np.sin((ori_360_m[i] / 180) * np.pi)),
               int(pos_y_m[i] + 60 * np.cos((ori_360_m[i] / 180) * np.pi)))
        desired_female_pos.append(list(ep1))

    dist_x_m = np.subtract(pos_x_m[2:], pos_x_m[:-2])
    dist_x_m = np.pad(dist_x_m, [1,1], 'constant')
    dist_y_m = np.subtract(pos_y_m[2:], pos_y_m[:-2])
    dist_y_m = np.pad(dist_y_m, [1,1], 'constant')
    dist_x_f = np.subtract(pos_x_f[2:], pos_x_f[:-2])
    dist_x_f = np.pad(dist_x_f, [1,1], 'constant')
    dist_y_f = np.subtract(pos_y_f[2:], pos_y_f[:-2])
    dist_y_f = np.pad(dist_y_f, [1,1], 'constant')

    dist_m = np.array([math.sqrt(x ** 2 + y ** 2) for x, y in zip(dist_x_m, dist_y_m)])
    speed_m = velocity(dist_m, timeperiod, timefactor)
    speed_m = np.multiply(speed_m, 60 / 1024)
    speed_m = np.median(np.lib.stride_tricks.sliding_window_view(speed_m, (window,)), axis=1)
    if window % 2 == 0:
        speed_m = np.pad(speed_m, [window//2 - 1, window//2], 'constant')
    else:
        speed_m = np.pad(speed_m, [window//2, window//2], 'constant')
    loco_dir_m = [(((np.arctan2(-y, x) / np.pi) * 180) + 90) % 360 for x, y in zip(dist_x_m, dist_y_m)]

    dist_f = np.array([math.sqrt(x ** 2 + y ** 2) for x, y in zip(dist_x_f, dist_y_f)])
    speed_f = velocity(dist_f, timeperiod, timefactor)
    speed_f = np.multiply(speed_f, 60 / 1024)
    speed_f = np.median(np.lib.stride_tricks.sliding_window_view(speed_f, (window,)), axis=1)
    if window % 2 == 0:
        speed_f = np.pad(speed_f, [window//2 - 1, window//2], 'constant')
    else:
        speed_f = np.pad(speed_f, [window//2, window//2], 'constant')
    loco_dir_f = [(((np.arctan2(-y, x) / np.pi) * 180) + 90) % 360 for x, y in zip(dist_x_f, dist_y_f)]

    return ori_smooth_m, ori_smooth_f, loco_dir_m, loco_dir_f, dist_m, dist_f, speed_m, speed_f, ang_vel_smooth_m, ang_vel_smooth_f


def compute_velocities_150(timeperiod, pos_x, pos_y, ori, no_of_trials, trial_period=None, timefactor=1e6, window=5):
    dist_x = np.subtract(pos_x[2:], pos_x[:-2])
    dist_x = np.pad(dist_x, [1,1], 'constant')
    dist_y = np.subtract(pos_y[2:], pos_y[:-2])
    dist_y = np.pad(dist_y, [1,1], 'constant')
    dist = np.array([math.sqrt(x ** 2 + y ** 2) for x, y in zip(dist_x, dist_y)])
    speed = velocity(dist, timeperiod, timefactor)

    speed = np.multiply(speed, 60 / 1024)
    speed_smooth = np.median(np.lib.stride_tricks.sliding_window_view(speed, (window,)), axis=1)

    if window % 2 == 0:
        speed_smooth = np.pad(speed_smooth, [window//2 - 1, window//2], 'constant')
    else:
        speed_smooth = np.pad(speed_smooth, [window//2, window//2], 'constant')

    speed_mean_smooth = np.convolve(speed_smooth, np.ones(window) / window, 'same')
    speed = np.array(speed_mean_smooth)

    
    turns = np.subtract((ori[2:]), (ori[:-2]))
    turns = np.pad(turns, [1,1], 'constant')
    turns_smooth = np.median(np.lib.stride_tricks.sliding_window_view(turns, (window,)), axis=1)
    if window % 2 == 0:
        turns_smooth = np.pad(turns_smooth, [window//2 - 1, window//2], 'constant')
    else:
        turns_smooth = np.pad(turns_smooth, [window//2, window//2], 'constant')
    turns_mean_smooth = np.convolve(turns_smooth, np.ones(5) / 5, 'same')
    ori_smooth = np.cumsum(turns_smooth) 

    ori_smooth = np.cumsum(turns_smooth)
    # ori_smooth = list(ori_smooth)
    ori_360_smooth = np.remainder(ori_smooth, 360)
    # angular_speed = velocity(list(turns_mean_smooth), timeperiod)

    if trial_period is None:
        pass
    else:
        for i in range(1, no_of_trials):
            turns[i * trial_period] = turns[(i * trial_period) - 1]

    angular_speed = velocity(turns, timeperiod, timefactor)
    angular_speed_smooth = np.median(np.lib.stride_tricks.sliding_window_view(angular_speed, (window,)), axis=1)
    if window % 2 == 0:
        angular_speed_smooth = np.pad(angular_speed_smooth, [window//2 - 1, window//2], 'constant')
    else:
        angular_speed_smooth = np.pad(angular_speed_smooth, [window//2, window//2], 'constant')
    angular_speed_mean_smooth = np.convolve(angular_speed_smooth, np.ones(5) / 5, 'same')

    angular_speed = np.array(angular_speed_mean_smooth)

    return ori_smooth, ori_360_smooth, angular_speed, speed


def compute_velocities(timeperiod, pos_x, pos_y, ori, no_of_trials, trial_period=None, timefactor=1e6):
    dist_x = np.diff(pos_x)
    dist_x = np.insert(dist_x, 0, 0)
    dist_y = np.diff(pos_y)
    dist_y = np.insert(dist_y, 0, 0)
    turns = np.diff(ori)
    turns = np.insert(turns, 0, 0)

    dist_x = list(dist_x)
    dist_y = list(dist_y)
    dist = []
    for i in range(len(dist_x)):
        dist.append(math.sqrt(dist_x[i] ** 2 + dist_y[i] ** 2))
    speed = velocity(dist, list(timeperiod), timefactor)
    speed = np.multiply(speed, 60 / 1024)
    speed_smooth = pd.Series(speed)
    speed_smooth = speed_smooth.rolling(3).median()
    speed_smooth[0:2] = [0, 0]
    speed_mean_smooth = np.convolve(speed_smooth, np.ones(3) / 3, 'valid')
    speed_mean_smooth = np.insert(speed_mean_smooth, 0, 0)
    speed_mean_smooth = np.append(speed_mean_smooth, 0)
    speed = np.array(speed_mean_smooth)

    turns = pd.Series(turns)
    turns_smooth = turns.rolling(3).median()
    turns_smooth[0:2] = [0, 0]
    turns_smooth = list(turns_smooth)
    turns_mean_smooth = np.convolve(turns_smooth, np.ones(3) / 3, 'valid')
    turns_mean_smooth = list(np.insert(turns_mean_smooth, 0, 0))
    turns_mean_smooth.append(0)

    ori_smooth = np.cumsum(turns_smooth)
    ori_smooth = np.insert(ori_smooth, 0, 0)
    ori_smooth = list(ori_smooth)
    ori_360_smooth = np.remainder(ori_smooth, 360)
    # angular_speed = velocity(list(turns_mean_smooth), timeperiod)

    if trial_period is None:
        pass
    else:
        for i in range(1, no_of_trials):
            turns[i * trial_period] = turns[(i * trial_period) - 1]

    angular_speed = velocity(list(turns), list(timeperiod), timefactor)
    angular_speed_smooth = pd.Series(angular_speed)
    angular_speed_smooth = angular_speed_smooth.rolling(3).median()
    angular_speed_smooth[0:2] = [0, 0]
    angular_speed_mean_smooth = np.convolve(angular_speed_smooth, np.ones(3) / 3, 'valid')
    angular_speed_mean_smooth = np.insert(angular_speed_mean_smooth, 0, 0)
    angular_speed_mean_smooth = np.append(angular_speed_mean_smooth, 0)

    angular_speed = np.array(angular_speed_mean_smooth)

    return ori_smooth, ori_360_smooth, angular_speed, speed


def compute_velocities_SPARC(timeperiod, pos_x, pos_y, ori, no_of_trials, trial_period):
    dist_x = np.diff(pos_x)
    dist_x = np.insert(dist_x, 0, 0)
    dist_y = np.diff(pos_y)
    dist_y = np.insert(dist_y, 0, 0)
    turns = np.diff(ori)
    turns = np.insert(turns, 0, 0)
    turns[turns > 20] = 0
    turns[turns < -20] = 0

    distance_to_px_ratio = 25/736

    dist_x = list(dist_x)
    dist_y = list(dist_y)
    dist = []
    for i in range(len(dist_x)):
        dist.append(math.sqrt(dist_x[i] ** 2 + dist_y[i] ** 2))
    speed = velocity(dist, list(timeperiod), 1e6)
    speed = np.multiply(speed, distance_to_px_ratio)
    speed_smooth = pd.Series(speed)
    speed_smooth = speed_smooth.rolling(3).median()
    speed_smooth[0:2] = [0, 0]
    speed_smooth = list(speed_smooth)
    speed_mean_smooth = np.convolve(speed_smooth, np.ones(3) / 3, 'valid')
    speed_mean_smooth = list(np.insert(speed_mean_smooth, 0, 0))
    speed_mean_smooth.append(0)

    # dist_x = pd.Series(dist_x)
    # dist_y = pd.Series(dist_y)
    # dist_x_smooth = dist_x.rolling(3).median()  # len(a) == len(a.rolling(n).median()
    # dist_x_smooth[0:2] = [0, 0]
    # dist_x_smooth = np.array(dist_x_smooth)
    # dist_y_smooth = dist_y.rolling(3).median()
    # dist_y_smooth[0:2] = [0, 0]
    # dist_y_smooth = np.array(dist_y_smooth)
    # dist_x_final = list(dist_x_smooth)
    # dist_y_final = list(dist_y_smooth)
    # dist = []
    # dist_smooth = []
    # dist_x = list(dist_x)
    # dist_y = list(dist_y)
    # for i in range(len(dist_x_smooth)):
    #     dist_smooth.append(math.sqrt(dist_x_smooth[i] ** 2 + dist_y_smooth[i] ** 2))
    # speed = velocity(dist_smooth, timeperiod)
    # speed = np.multiply(speed, 60 / 1024)
    speed = np.array(speed_mean_smooth)

    turns = pd.Series(turns)
    turns_smooth = turns.rolling(3).median()
    turns_smooth[0:2] = [0, 0]
    turns_smooth = list(turns_smooth)
    turns_mean_smooth = np.convolve(turns_smooth, np.ones(3) / 3, 'valid')
    turns_mean_smooth = list(np.insert(turns_mean_smooth, 0, 0))
    turns_mean_smooth.append(0)

    ori_smooth = np.cumsum(turns_smooth)
    ori_smooth = np.insert(ori_smooth, 0, 0)
    ori_smooth = list(ori_smooth)
    ori_360_smooth = np.remainder(ori_smooth, 360)
    # angular_speed = velocity(list(turns_mean_smooth), timeperiod)

    for i in range(1, no_of_trials):
        turns[i * trial_period] = turns[(i * trial_period) - 1]
    angular_speed = velocity(list(turns), list(timeperiod), 1e6)
    angular_speed_smooth = pd.Series(angular_speed)
    angular_speed_smooth = angular_speed_smooth.rolling(3).median()
    angular_speed_smooth[0:2] = [0, 0]
    angular_speed_smooth = list(angular_speed_smooth)
    angular_speed_mean_smooth = np.convolve(angular_speed_smooth, np.ones(3) / 3, 'valid')
    angular_speed_mean_smooth = list(np.insert(angular_speed_mean_smooth, 0, 0))
    angular_speed_mean_smooth.append(0)

    angular_speed = np.array(angular_speed_mean_smooth)

    speed_change = np.diff(angular_speed_mean_smooth)
    speed_change = np.insert(speed_change, 0, 0)
    angular_acceleration = velocity(list(speed_change), list(timeperiod), 1e6)

    return ori_smooth, ori_360_smooth, angular_speed, speed, angular_acceleration

def locomotion_direction(pos_x, pos_y, correct_ori, window=3):
    """
    Get the direction of displacement of the fly for each frame
    :param pos_x: x position of the fly
    :param pos_y: y position of the fly
    :param correct_ori: corrected orientation of the fly, range -180, 180
    :param window: window for smoothing
    :return: direction of displacement of the fly for each frame, range 0, 360
    :rtype: np.array
    :return: distance of the fly for each frame
    :rtype: np.array
    """
    if pos_x.shape[0] % window != 0 :
        pos_x = pos_x[:(pos_x.shape[0]//window) * window]
        pos_y = pos_y[:(pos_y.shape[0] // window) * window]

    if window%2==0:
        dist_x = np.pad(np.subtract(pos_x[window:], pos_x[:-window]), (window, 0), mode='constant', constant_values=0)
        dist_y = np.pad(np.subtract(pos_y[window:], pos_y[:-window]), (window, 0), mode='constant', constant_values=0)
    else:
        dist_x = np.pad(np.subtract(pos_x[window:], pos_x[:-window]), (window, 0), mode='constant', constant_values=0)
        dist_y = np.pad(np.subtract(pos_y[window:], pos_y[:-window]), (window, 0), mode='constant', constant_values=0)
    dist = np.sqrt(np.add(dist_x**2, dist_y**2))
    dir = (((np.arctan2(dist_y, dist_x)/np.pi) * 180) - 90)%360
    ##
    indices = np.argwhere(dist<2).flatten()
    dir[indices] = correct_ori[indices]
    return dir, dist

def speed_components(speed, heading_direction, ori_list):
    """
    Get the parallel and perpendicular components of the speed w.r.t the heading direction
    :param speed: speed of the fly
    :type speed: np.array
    :param heading_direction: heading direction of the fly
    :type heading_direction: np.array
    :param ori_list: corrected orientation of the fly
    :type ori_list: np.array
    :return: parallel and perpendicular components of the speed w.r.t the heading direction
    :rtype: (np.array, np.array)
    """
    ## to make sure that the direction is right, the head-tail direction has to be specified
    relative_heading_direction = (((heading_direction - ori_list)%360)/180)*np.pi
    v_parallel = speed * np.cos(relative_heading_direction)
    v_perpendicular = speed * np.sin(relative_heading_direction)
    return v_parallel, v_perpendicular