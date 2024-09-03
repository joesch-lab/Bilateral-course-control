import numpy as np


def remove_bad_trials(data, p=2, strict=False, non_stim=False):
    new_data = []
    ## data[0] should always be opto_response
    n = data[0].shape[1]
    trial_quality = data[0][:, p+3]
    if not strict:
        if not non_stim:
            for i in range(len(data)):
                if data[i].ndim == 2:
                    indices = np.where(trial_quality != 1)
                    new_data.append(data[i][indices])
                elif data[i].ndim == 3:
                    indices = np.where(trial_quality != 1)
                    new_data.append(data[i][indices])
        else:
            for i in range(len(data)):
                if data[i].ndim == 2:
                    indices = np.where((trial_quality != -1) & (trial_quality != 1))
                    new_data.append(data[i][indices])
                elif data[i].ndim == 3:
                    indices = np.where((trial_quality != -1) & (trial_quality != 1))
                    new_data.append(data[i][indices])
    if strict:
        if not non_stim:
            for i in range(len(data)):
                if data[i].ndim == 2:
                    indices = np.where((trial_quality != -1) & (trial_quality != 1) & (trial_quality != 5))
                    new_data.append(data[i][indices])
                elif data[i].ndim == 3:
                    indices = np.where((trial_quality != -1) & (trial_quality != 1) & (trial_quality != 5))
                    new_data.append(data[i][indices])
        else:
            for i in range(len(data)):
                if data[i].ndim == 2:
                    indices = np.where((trial_quality != -1) & (trial_quality != 1) & (trial_quality != 5) & (trial_quality != -5))
                    new_data.append(data[i][indices])
                elif data[i].ndim == 3:
                    indices = np.where((trial_quality != -1) & (trial_quality != 1) & (trial_quality != 5) & (trial_quality != -5))
                    new_data.append(data[i][indices])
    return new_data
