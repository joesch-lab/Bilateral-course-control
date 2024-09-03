import numpy as np

def reject_nonresponsive_trials(opto_response_angular_speed, opto_response_speed, p=3, stim_period=600):
    """
    remove trials in which the fly does not move at all based on the angular speed and speed data
    :param opto_response_angular_speed: The angular speed data.
    :type opto_response_angular_speed: numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param opto_response_speed: The speed data.
    :type opto_response_speed: numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param p: The number of independent sections in the stimulus, required for parsing the np array, can be taken out.
    :type p: int
    :param stim_period: The period after which the rotation starts, in frames, 5s at 60Hz is 300 frames.
    :type stim_period: int
    :returns: The indices of the trials that are responsive (new_indices) and non-responsive (rejected_indices).
    :rtype: (numpy.ndarray, numpy.ndarray)
    """
    speed_mean = np.nanmean(opto_response_speed[:, p+4:], axis=1)
    angular_speed_var = np.std(opto_response_angular_speed[:, p+4+stim_period:], axis=1)
    rejected_indices = np.intersect1d(np.argwhere(angular_speed_var < 50).flatten(), np.argwhere(speed_mean < 1).flatten())
    # angular_speed_var = np.nanmean(np.abs(opto_response_angular_speed[:, p + 4 + stim_period:]), axis=1)
    # rejected_indices = np.union1d(np.argwhere(angular_speed_var < 5).flatten(), np.argwhere(speed_mean < 1).flatten())
    original_indices = np.arange(0, opto_response_angular_speed.shape[0])
    new_indices = np.delete(original_indices, rejected_indices)
    return new_indices, rejected_indices
