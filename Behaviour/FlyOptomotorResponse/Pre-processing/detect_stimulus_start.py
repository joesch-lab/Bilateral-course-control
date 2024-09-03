import numpy as np

def stimulus_pattern(my_list, check = 0):
    """
    returns the frame number where the stimulus changes based on the default value of the stimulus
    :param my_list: list of stimulus values
    :type my_list: list
    :param check: value of stimulus before a change, the default value of my_list
    :type check: int 
    :return: list of frame numbers where the stimulus changes and the new stimulus value
    :rtype: (list, list)
    """
    frame_changes = []
    stim_change = []
    for i in range(2,len(my_list)):
        if my_list[i-1] == check and my_list[i-2] == check and my_list[i] != check:
            frame_changes.append(i)
            stim_change.append(my_list[i])
    return frame_changes, stim_change


def stimulus_pattern_hemi(direction_left, direction_right, direction_center, direction_rear, check = 0):
    """
    returns the frame number where the stimulus changes for the hemifield stimulus
    :param direction_left, direction_right, direction_center, direction_rear: arrays of stimulus values
    :type direction_left, direction_right, direction_center, direction_rear: np.ndarray
    :param check: value of stimulus that is considered as change
    :type check: int 
    :return: list of frame numbers where the stimulus changes
    :rtype: list
    """
    frame_changes = []
    stim_change = []
    my_list = np.stack((direction_left, direction_right, direction_center, direction_rear))

    for i in range(2,my_list.shape[1]):
        if np.all(my_list[:, i-1] == np.ones(4)*check) and (not(np.all(my_list[:, i] == np.ones(4)*check))):
            frame_changes.append(i)
            # stim_change.append(my_list[i])
    return frame_changes


def stimulus_pattern_specific(my_list, check):
    """
    returns the frame number where the stimulus changes based on what stimulus value is considered as change
    :param my_list: list of stimulus values
    :type my_list: list
    :param check: value of stimulus that is considered as change
    :type check: int 
    :return: list of frame numbers where the stimulus changes and the new stimulus value
    :rtype: (list, list)
    """
    frame_changes = []
    stim_change = []
    for i in range(1,len(my_list)):
        if my_list[i-1] != my_list[i] and my_list[i] == check:
            frame_changes.append(i)
            stim_change.append(my_list[i])
    return frame_changes,stim_change


def stimulus_pattern_new(my_list):
    """
    returns the frame number where the stimulus changes
    :param my_list: list of stimulus values
    :type my_list: list
    :return: list of frame numbers where the stimulus changes and the new stimulus value
    :rtype: (list, list)
    """
    frame_changes = []
    stim_change = []
    frame_changes.append(0)
    stim_change.append(my_list[0])
    for i in range(1,len(my_list)):
        if my_list[i-1] != my_list[i]:
            frame_changes.append(i)
            stim_change.append(my_list[i])
    return frame_changes,stim_change


def stimulus_pattern_change(my_list):
    """
    returns the frame number where the stimulus changes but the change is from 0 to a non-zero value
    :param my_list: list of stimulus values
    :type my_list: list
    :return: list of frame numbers where the stimulus changes and the new stimulus value
    :rtype: (list, list)
    """    
    frame_changes = []
    stim_change = []
    if my_list[0] != 0:
        frame_changes.append(0)
        stim_change.append(my_list[0])
    for i in range(1, len(my_list)):
        if my_list[i] != 0:
            if my_list[i-1] != my_list[i]:
                frame_changes.append(i)
                stim_change.append(my_list[i])
    return frame_changes, stim_change


def stimulus_pattern_off(my_list, check=0):
    """
    returns the frame number where the stimulus changes based on when stimulus value is not 0 and equal to check
    :param my_list: list of stimulus values
    :type my_list: list
    :param check: value of stimulus that is considered as change
    :type check: int 
    :return: list of frame numbers where the stimulus changes and the new stimulus value
    :rtype: list
    """
    frame_changes = []
    for i in range(1, len(my_list)):
        if my_list[i-1] != 0 and my_list[i] == check:
            frame_changes.append(i)
    return frame_changes