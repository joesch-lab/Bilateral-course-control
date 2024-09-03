import copy
import numpy as np


def remove_saccades(opto_response_angular_speed, opto_response, saccade_peaks, window=[5, 5], p=2, 
                    remove='anti', stim_period=600, filler_type='zero'):
    """
    remove saccades from the cumulative turn and angular speed
    output to be used to compare with origina data to see the effect of smooth tracking
    """
    new_opto_response_angular_speed = copy.deepcopy(opto_response_angular_speed) ## copy because we are going to change this array
                                                                                ## and we want to preserve the original array
    opto_turns = np.apply_along_axis(lambda x: np.diff(x), axis=1, arr=opto_response[:, p+4:])
    opto_turns = np.insert(opto_turns, 0, opto_response[:, p+4], 1)
    new_opto_response_angular_speed = new_opto_response_angular_speed[:, p+3:]

    for i in range(opto_turns.shape[0]):
        if remove == 'anti' or remove == 'all':
            if saccade_peaks[i, 0] == -1:
                anti_indices = np.argwhere(saccade_peaks[i, p+1+stim_period:] < 0)
                for index in anti_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = np.zeros(window[0] + window[1])
                        if filler_type == 'zero':
                            filler = np.zeros(window[0] + window[1])
                        elif filler_type == 'linear':
                            filler = np.interp(np.arange(0, window[0] + window[1], 1), [-1, window[0] + window[1] + 1],
                                               [new_opto_response_angular_speed[i, index[0] + stim_period - window[0] - 1],
                                                new_opto_response_angular_speed[i, index[0] + stim_period - window[0] + 1]])
                            filler = np.random.rand(window[0] + window[1]) * 25 + filler
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = filler
            if saccade_peaks[i, 0] == 1:
                anti_indices = np.argwhere(saccade_peaks[i, p+1+stim_period:] > 0)
                for index in anti_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = np.zeros(window[0] + window[1])
                        if filler_type == 'zero':
                            filler = np.zeros(window[0] + window[1])
                        elif filler_type == 'linear':
                            filler = np.interp(np.arange(0, window[0] + window[1], 1), [-1, window[0] + window[1] + 1],
                                               [new_opto_response_angular_speed[i, index[0] + stim_period - window[0] - 1],
                                                new_opto_response_angular_speed[i, index[0] + stim_period - window[0] + 1]])
                            filler = np.random.rand(window[0] + window[1]) * 25 + filler
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = filler
        if remove == 'syn' or remove == 'all':
            if saccade_peaks[i, 0] == 1:
                syn_indices = np.argwhere(saccade_peaks[i, p+1 + stim_period:] < 0)
                for index in syn_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        opto_turns[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]] = np.zeros(window[0] + window[1])
                        if filler_type == 'zero':
                            filler = np.zeros(window[0] + window[1])
                        elif filler_type == 'linear':
                            filler = np.interp(np.arange(0, window[0] + window[1], 1), [-1, window[0] + window[1] + 1],
                                               [new_opto_response_angular_speed[i, index[0] + stim_period - window[0] - 1],
                                                new_opto_response_angular_speed[i, index[0] + stim_period - window[0] + 1]])
                            filler = np.random.rand(window[0] + window[1]) * 25 + filler
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = filler
            if saccade_peaks[i, 0] == -1:
                syn_indices = np.argwhere(saccade_peaks[i, p+1 + stim_period:] > 0)
                for index in syn_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        opto_turns[i, index[0] + stim_period-window[0]:index[0] + stim_period + window[1]] = np.zeros(window[0] + window[1])
                        if filler_type == 'zero':
                            filler = np.zeros(window[0] + window[1])
                        elif filler_type == 'linear':
                            filler = np.interp(np.arange(0, window[0] + window[1], 1), [-1, window[0] + window[1] + 1],
                                               [new_opto_response_angular_speed[i, index[0] + stim_period - window[0] - 1],
                                                new_opto_response_angular_speed[i, index[0] + stim_period - window[0] + 1]])
                            filler = np.random.rand(window[0] + window[1]) * 25 + filler
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = filler
    new_optomotor_response = np.apply_along_axis(lambda x: np.cumsum(x), axis=1, arr=opto_turns)
    for i in range(p+4):
        new_optomotor_response = np.insert(new_optomotor_response, i, opto_response[:, i], 1)
    for i in range(p+3):
        new_opto_response_angular_speed = np.insert(new_opto_response_angular_speed, i, new_opto_response_angular_speed[:, i], 1)
    return new_optomotor_response, new_opto_response_angular_speed


def keep_saccades(opto_response_angular_speed, opto_response, saccade_peaks, window=[5, 5], p=2, keep='anti', stim_period=600):
    """
    keep only saccadic turns and remove smooth turning in the cumulative turn and angular speed data
    """
    new_opto_response_angular_speed = np.zeros(opto_response_angular_speed.shape)
    opto_turns = np.apply_along_axis(lambda x: np.diff(x), axis=1, arr=opto_response[:, p+4:])
    opto_turns = np.insert(opto_turns, 0, opto_response[:, p+4], 1)
    new_opto_turns = np.zeros(opto_turns.shape)
    new_opto_response_angular_speed = new_opto_response_angular_speed[:, p+3:]

    for i in range(opto_turns.shape[0]):
        if keep == 'anti' or keep == 'all':
            if saccade_peaks[i, 0] == -1:
                anti_indices = np.argwhere(saccade_peaks[i, p+1+stim_period:] < 0)
                for index in anti_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        new_opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]]
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = \
                            opto_response_angular_speed[i, p+3+index[0] + stim_period - window[0]:p+3+index[0] + stim_period + window[1]]
            if saccade_peaks[i, 0] == 1:
                anti_indices = np.argwhere(saccade_peaks[i, p+1+stim_period:] > 0)
                for index in anti_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        new_opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]]
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = \
                            opto_response_angular_speed[i, p+3+index[0] + stim_period - window[0]:p+3+index[0] + stim_period + window[1]]
        if keep == 'syn' or keep == 'all':
            if saccade_peaks[i, 0] == 1:
                syn_indices = np.argwhere(saccade_peaks[i, p+1 + stim_period:] < 0)
                for index in syn_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        new_opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]]
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = \
                            opto_response_angular_speed[i, p+3+index[0] + stim_period - window[0]:p+3+index[0] + stim_period + window[1]]
            if saccade_peaks[i, 0] == -1:
                syn_indices = np.argwhere(saccade_peaks[i, p+1 + stim_period:] > 0)
                for index in syn_indices:
                    if index[0] > window[0] and index[0] < (opto_turns.shape[1]-stim_period) - window[1]:
                        new_opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = opto_turns[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]]
                        new_opto_response_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = \
                            opto_response_angular_speed[i, p+3+index[0] + stim_period - window[0]:p+3+index[0] + stim_period + window[1]]
    new_optomotor_response = np.apply_along_axis(lambda x: np.cumsum(x), axis=1, arr=new_opto_turns)
    for i in range(p+4):
        new_optomotor_response = np.insert(new_optomotor_response, i, opto_response[:, i], 1)
    for i in range(p+3):
        new_opto_response_angular_speed = np.insert(new_opto_response_angular_speed, i, new_opto_response_angular_speed[:, i], 1)
    return new_optomotor_response, new_opto_response_angular_speed