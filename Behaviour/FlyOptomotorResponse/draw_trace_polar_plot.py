import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter

import matplotlib.pyplot as plt

import matplotlib
plt.style.use('roshan_style.mplstyle')
figsize = matplotlib.rcParams['figure.figsize']

def draw_trace_polar_plot(data_list, Dir, filename, labels, fly_indices=None):
    """
    Draw a polar plot for the provided data, used for plotting fly orientation to visualise turning response. Bad old code.
    :param data_list: The data to be plotted
    :type data_list: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    """
    n, m = subplot_arrangement(len(data_list))
    fig, ax = plt.subplots(n, m, figsize=(figsize[0]* m, figsize[1]* n), subplot_kw={'projection':"polar"})
    fig2, ax2 = plt.subplots(n, m, figsize=(figsize[0] * m, figsize[1] * n), subplot_kw={'projection': "polar"})
    fig3, ax3 = plt.subplots(n, m, figsize=(figsize[0] * m, figsize[1] * n), subplot_kw={'projection': "polar"})
    end = 60 + 75
    for i in range(len(data_list)):
        if ('FlpD' in filename) and (i == 0):
            data_list_smoothed = savgol_filter(gaussian_filter1d(data_list[i][:, 240:240+end+1], sigma=3, axis=1), 31, 2)
        else:
            data_list_smoothed = gaussian_filter1d(data_list[i][:, 240:240 + end + 1], sigma=3, axis=1)
        data_list_smoothed = ((data_list_smoothed/180) * np.pi) + np.pi/2
        r = np.arange(0, end+1)

        ax[i].plot(np.mean(data_list_smoothed, axis=0), r, lw=2, alpha=0.9, c='k')
        ax2[i].plot(np.mean(data_list_smoothed, axis=0), r, lw=2, alpha=0.9, c='k')
        print(data_list_smoothed.shape)

        if fly_indices != None:
            for ids in np.unique(fly_indices[i]):
                indices = np.argwhere(fly_indices[i] == ids).flatten()
                if indices.shape[0] < 5:
                    pass
                else:
                    # print(ids, indices.shape)
                    ax2[i].plot(np.nanmean(data_list_smoothed[indices], axis=0), r, lw=0.5, alpha=0.5, c='grey')
                    for j in range(indices.shape[0]):
                        ax3[i].plot(data_list_smoothed[indices][j], r, lw=0.4, alpha=0.3, c='grey')

        max_id = np.random.choice(np.arange(10))
        indices = np.argwhere(fly_indices[i] == max_id).flatten()
        for j in range(min(data_list_smoothed[indices].shape[0], 10)):
            ax[i].plot(data_list_smoothed[indices][j], r, lw=0.4, alpha=0.3, c='grey')

        ax[i].set_title(labels[i])
        ax[i].grid(True)
        ax[i].set_rticks([0, 60, end])
        ax[i].set_yticklabels([])
        ax[i].set_xticks([])
        ax[i].spines['polar'].set_visible(False)
        ax[i].fill_between(np.arange(0, 2, 1./180)*np.pi, 0, 60, alpha=0.1, color='grey')

        ax2[i].set_title(labels[i])
        ax2[i].grid(True)
        ax2[i].set_rticks([0, 60, end])
        ax2[i].set_yticklabels([])
        ax2[i].set_xticks([])
        ax2[i].spines['polar'].set_visible(False)
        ax2[i].fill_between(np.arange(0, 2, 1./180)*np.pi, 0, 60, alpha=0.1, color='grey')

        ax3[i].set_title(labels[i])
        ax3[i].grid(True)
        ax3[i].set_rticks([0, 60, end])
        ax3[i].set_yticklabels([])
        ax3[i].set_xticks([])
        ax3[i].spines['polar'].set_visible(False)
        ax3[i].fill_between(np.arange(0, 2, 1./180)*np.pi, 0, 60, alpha=0.1, color='grey')

    fig.suptitle(filename.split('_')[0], fontsize=16)
    fig2.suptitle(filename.split('_')[0], fontsize=16)
    fig3.suptitle(filename.split('_')[0], fontsize=16)
    # fig.savefig(os.path.join(Dir, filename+'polarplot.pdf'), dpi=600, transparent=True)
    # fig.savefig(os.path.join(Dir, filename + 'polarplot.svg'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir, filename + 'polarplot.png'))
    fig2.savefig(os.path.join(Dir, filename + 'polarplot2.png'))
    fig3.savefig(os.path.join(Dir, filename + 'polarplot3.png'))
    return 1