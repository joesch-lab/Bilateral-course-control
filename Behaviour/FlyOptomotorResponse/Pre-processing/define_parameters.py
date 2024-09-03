import numpy as np

def define_parameters(stim_left, stim_right, stim_center, contrast_left, contrast_right, stimulus, freq_left, freq_right, PARAMETERS=['frequency', 'bars'], new=True):
    if new:
        n=1
    else:
        n=2
    parameter=[]
    for i in range(len(PARAMETERS)):
        if PARAMETERS[i] == 'stim_left':
            parameter.append(np.array(stim_left))
        elif PARAMETERS[i] == 'frequency_left':
            parameter.append(np.repeat(freq_left, n))
        elif PARAMETERS[i] == 'frequency_right':
            parameter.append(np.repeat(freq_right, n))
        elif PARAMETERS[i] == 'contrast_left':
            parameter.append(np.repeat(contrast_left, n))
        elif PARAMETERS[i] == 'contrast_right':
            parameter.append(np.repeat(contrast_right, n))
        elif PARAMETERS[i] == 'stim_right':
            parameter.append(np.array(stim_right))
        elif PARAMETERS[i] == 'stim_center':
            parameter.append(np.array(stim_center))
        elif PARAMETERS[i] == 'stimulus' or PARAMETERS[i] == 'density':
            parameter.append(np.repeat(stimulus, n))
            # parameter.append(stimulus)
    return parameter


def define_parameters2(stim_left, stim_right, contrast_left, contrast_right, stimulus, freq_left, freq_right, PARAMETERS=['frequency', 'bars']):
    parameter=[]
    for i in range(len(PARAMETERS)):
        if PARAMETERS[i] == 'stim_left':
            parameter.append(np.array(stim_left))
        elif PARAMETERS[i] == 'frequency_left':
            parameter.append(np.repeat(freq_left, 2))
        elif PARAMETERS[i] == 'frequency_right':
            parameter.append(np.repeat(freq_right, 2))
        elif PARAMETERS[i] == 'contrast_left':
            parameter.append(np.repeat(contrast_left, 2))
        elif PARAMETERS[i] == 'contrast_right':
            parameter.append(np.repeat(contrast_right, 2))
        elif PARAMETERS[i] == 'stim_right':
            parameter.append(np.array(stim_right))
        elif PARAMETERS[i] == 'stimulus' or PARAMETERS[i] == 'density':
            parameter.append(np.repeat(stimulus, 2))
            # parameter.append(stimulus)
    return parameter


def define_parameters4(stim_left, stim_right, stim_center, stim_rear, contrast_left, contrast_right, stimulus, freq_left, freq_right, PARAMETERS=['frequency', 'bars']):
    parameter=[]
    for i in range(len(PARAMETERS)):
        if PARAMETERS[i] == 'stim_left':
            parameter.append(np.array(stim_left))
        elif PARAMETERS[i] == 'frequency_left':
            parameter.append(np.repeat(freq_left, 2))
        elif PARAMETERS[i] == 'frequency_right':
            parameter.append(np.repeat(freq_right, 2))
        elif PARAMETERS[i] == 'contrast_left':
            parameter.append(np.repeat(contrast_left, 2))
        elif PARAMETERS[i] == 'contrast_right':
            parameter.append(np.repeat(contrast_right, 2))
        elif PARAMETERS[i] == 'stim_right':
            parameter.append(np.array(stim_right))
        elif PARAMETERS[i] == 'stimulus' or PARAMETERS[i] == 'density':
            parameter.append(np.repeat(stimulus, 2))
        elif PARAMETERS[i] == 'stim_center':
            parameter.append(np.array(stim_center))
        elif PARAMETERS[i] == 'stim_rear':
            parameter.append(np.array(stim_rear))
            # parameter.append(stimulus)
    return parameter
