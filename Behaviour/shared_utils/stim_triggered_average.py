import numpy as np

def stim_triggered_average(frame_changes, data, pre_stim_time= 60, post_stim_time=60):
    """
    split data into trials based on beginning of stimulus given by frame_changes
    pre and post-stimulus period given by the respective values
    :param frame_changes: list of frame numbers where the stimulus changes
    :type frame_changes: list or np.ndarray
    :param data: list of data values
    :type data: list or np.ndarray
    :param pre_stim_time: time before the stimulus starts
    :type pre_stim_time: int
    :param post_stim_time: time after the stimulus starts
    :type post_stim_time: int
    :return: data_list - array of data values for each trial of shape (n, m) where m=[pre_stim_time + post_stim_time + 1]
    :rtype: np.ndarray
    """
    try:
        if len(frame_changes)==0:
            return np.array([data])
    except:
        if frame_changes.shape==0:
            return np.array([data])
    data_list = []
    length = pre_stim_time + post_stim_time + 1
    for i in frame_changes:
        print(i-pre_stim_time, i+post_stim_time+1)
        if len(data[i-pre_stim_time:i+post_stim_time+1]) == length:
            data_list.append(np.array(data[i-pre_stim_time:i+post_stim_time+1]))
        else:
            data_list.append(np.zeros(pre_stim_time + post_stim_time + 1))
    data_list = np.array(data_list)

    return data_list
