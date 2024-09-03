import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import sem
from matplotlib import cm

import matplotlib.pyplot as plt

import matplotlib
plt.style.use('roshan_style.mplstyle')
figsize = matplotlib.rcParams['figure.figsize']

def draw_trace_each_fly(data_list, fly_indices, Dir, filename, labels, top=100, bottom=-100, cmap='tab20', smoothed='N', ls=[], turn_or_angular_speed=True,
                        axvline=120, start=180, end=-1):
    """
    Draw line plot showing the mean value for each fly for each stimulus condition.
    :param data_list: The data to be plotted
    :type data_list: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param fly_indices: The indices of the flies corresponding to each trial
    :type fly_indices: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n,) where n is the number of trials
    """
    ## data_list is a list of length n where n is the number of stimulus conditions
    ## fly_indices is a list length no arrays, where each array is indices for that stimulus condition
    n,m = subplot_arrangement(len(data_list))
    if not ls:
        ls = ['solid'] * len(labels)
    fig, ax = plt.subplots(n, m, figsize=(figsize[0]* m, figsize[1]* n)) ## plots mean trace of each fly
    fig1, ax1 = plt.subplots(n, m, figsize=(figsize[0]* m, figsize[1]* n)) ## plots the mean and sem for each fly
    # fig2, ax2 = plt.subplots(1, 5, figsize=(figsize[0]* 5, figsize[1]* 1)) ## plots the average traces of a few selected flies

    colors = [cm.get_cmap(cmap)(i) for i in np.linspace(0, 1, np.unique(fly_indices[0]).shape[0])]

    for i in range(len(data_list)):
        if smoothed=='Y':
            data_list_smoothed = gaussian_filter1d(data_list[i], sigma=3, axis=1)
        else:
            data_list_smoothed = np.array(data_list[i])

        if start == 0 and end == -1:
            pass
        else:
            print(data_list_smoothed.shape)
            data_list_smoothed = data_list_smoothed[:, start:end]

        l = data_list_smoothed.shape[1]
        if len(data_list) == 1:
            axis = ax
            axis1 = ax1
        else:
            if n == 1 or m == 1:
                axis = ax[i]
                axis1 = ax1[i]
            else:
                axis = ax[i // m, i % m]
                axis1 = ax1[i // m, i % m]
        response_values = []
        f = 0
        for k, ids in enumerate(np.unique(fly_indices[i])):
            indices = np.argwhere(fly_indices[i] == ids).flatten()
            axis.plot(np.nanmean(data_list_smoothed[indices], axis=0), lw=1, c='grey', label=str(ids))
            if turn_or_angular_speed:
                response_values.append(data_list_smoothed[indices][:, -1])
                axis1.text(ids - 0.6, np.nanmean(data_list_smoothed[indices][:, -1]), str(ids))
                axis1.scatter([ids] * data_list_smoothed[indices][:, -1].shape[0], data_list_smoothed[indices][:, -1], s=20, alpha=0.2, c='k')
            else:
                response_values.append(np.nanmean(data_list_smoothed[indices][:, stim_period:], axis=1))
                axis1.text(ids - 0.6, np.nanmean(np.nanmean(data_list_smoothed[indices][:, stim_period:], axis=1)), str(ids))
                axis1.scatter([ids] * data_list_smoothed[indices][:, -1].shape[0], np.nanmean(data_list_smoothed[indices][:, stim_period:], axis=1), s=20, alpha=0.2, c='k')
        axis1.errorbar(np.unique(fly_indices[i]), [np.nanmean(x) for x in response_values], yerr=[sem(x) for x in response_values], fmt='none', color='k', ecolor='k',elinewidth=2.5, capsize=0)
        axis1.scatter(np.unique(fly_indices[i]), [np.nanmean(x) for x in response_values], s=30, alpha=0.5, c='k')
        axis.plot(np.nanmean(data_list_smoothed, axis=0), lw=3, c='k')
        axis.set_ylim(top=top, bottom=bottom)
        axis.set_title(labels[i])
        axis.axhline(y=0, ls='--', lw=1, c='grey')
        axis.legend(loc='upper left', ncol=3)
        if i != 0:
            axis.set_yticks([])
            axis1.set_yticks([])
        axis1.set_xlim(left=np.amin(np.unique(fly_indices[i]))-1, right=np.amax(np.unique(fly_indices[i]))+1)
        axis1.set_xticks([])
        axis1.set_title(labels[i])
        axis1.set_ylim(top=top+200, bottom=bottom-200)
        axis1.axhline(y=0, lw=1, ls='--', c='gray')

        axis.set_ylim(top=top, bottom=bottom)
        axis.set(aspect=(data_list_smoothed.shape[1]) / (top - bottom))
        axis.axvline(x=axvline, lw=1, ls='--', c='grey')
        axis.axhline(y=0, lw=1, ls='--', c='grey', alpha=0.5)
        axis.set_ylim(top=top + 20, bottom=bottom)
        axis.set_yticks([top, 0, bottom])
        axis.set_xticks([0, axvline, data_list_smoothed.shape[1]])
        axis.set_xticklabels([0, int(axvline // 60), int(data_list_smoothed.shape[1] // 60)])
        axis.spines['bottom'].set_bounds(high=data_list_smoothed.shape[1], low=0)
        axis.spines['left'].set_bounds(high=top, low=bottom)
        axis.spines['left'].set_position(('outward', 4))
        axis.spines['bottom'].set_position(('outward', 2))
    i += 1
    while (i < n * m):
        fig.delaxes(ax[i // m, i % m])
        fig1.delaxes(ax1[i // m, i % m])
        i += 1
    # ax[i // m, i % m].legend()
    fig.savefig(r'{}/{}'.format(Dir,  filename+'fly_id.png'))
    fig1.savefig(r'{}/{}'.format(Dir, filename + 'fly_id_full_data.png'))
    fig.savefig(r'{}/{}'.format(Dir,  filename+'fly_id.svg'))
    fig1.savefig(r'{}/{}'.format(Dir, filename + 'fly_id_full_data.svg'))
    return 1