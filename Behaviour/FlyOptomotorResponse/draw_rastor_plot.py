import numpy as np
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import minmax_scale
from scipy.ndimage import gaussian_filter1d
from numpy.random import default_rng
from scipy.signal import savgol_filter
from scipy.stats import sem

import matplotlib.pyplot as plt

import matplotlib
plt.style.use('roshan_style.mplstyle')
figsize = matplotlib.rcParams['figure.figsize']

def draw_rastor_plot(data_list, Dir, filename, cmap=ListedColormap(all_color_data['PiYG7']), vmin=0, vmax=100, labels=[], smoothed='Y', normalize='N', sort=[0, 1], to_plot=0,
                     PDF='N', axvline=None, trace_average='N', start=0, end=-1, fly_indices=None):
    """
    Draw a rastor plot of the provided data, used for plotting speed or angular speed. Bad old code.
    :param data_list: The data to be plotted
    :type data_list: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :trace_average: 'Y' or 'N' if 'Y' then another axis is added at the bottom of the rastor and the mean trace is plotted
    :type trace_average: str
    :returns: mean and median of data and error (standard error of mean of the data), relevant for angular speed data only)
    :rtype: (numpy.ndarray, numpy.ndarray, numpy.ndarray)
    """
    if axvline is None:
        axvline=data_list[0].shape[1]//2
    n, m = subplot_arrangement(len(data_list)) # n rows and m columns
    if trace_average == 'N':
        fig, ax = plt.subplots(n, m, figsize=(figsize[0]* m, figsize[1]* n), gridspec_kw={'height_ratios': [1] * n, 'width_ratios': [1] * m})
    elif trace_average == 'Y':
        fig, ax = plt.subplots(n, m, figsize=(figsize[0] * m + 1, figsize[1] * n +1), gridspec_kw={'height_ratios': [1] * n, 'width_ratios': [1] * m, 'hspace':0.1, 'bottom':0.0})
    trace_ylim = [-20, 20]
    mean_trace = []
    mean_trace = []
    error_trace = []
    median_trace = []
    for i in range(len(data_list)):
        if smoothed=='Y':
            # data_list_smoothed = gaussian_filter1d(data_list[i], sigma=3, axis=1)
            data_list_smoothed = np.array(data_list[i])
        else:
            data_list_smoothed = np.array(data_list[i])

        if normalize=='Y':
            data_list_smoothed = minmax_scale(data_list_smoothed, (-1, 1), axis = 1)
            vmax = 1
            vmin = -1
        else:
            pass

        if len(sort)==2:
            indices = np.argsort(np.nanmean(data_list_smoothed[:, sort[0]:sort[1]], axis=1)).flatten()
            data_list_smoothed = data_list_smoothed[indices]
        else:
            if sort == [0]:
                ## new way of sorting
                average_response_rounded = []
                for start_i in np.arange(start, data_list_smoothed.shape[1], 30):
                    ## need to round the values to some integer to get lexsort to work nicely since you need same values
                    average_response_rounded.append(np.nanmean(data_list_smoothed[:, start_i: start_i + 30], axis=1)//50)
                indices = np.lexsort(average_response_rounded).flatten()
                data_list_smoothed = data_list_smoothed[indices]
                if fly_indices != None:
                    fly_indices[i] = fly_indices[i][indices]
            else:
                pass

        if start==0 and end==-1:
            pass
        else:
            data_list_smoothed = data_list_smoothed[:, start:end]

        l = data_list_smoothed.shape[1]
        if len(data_list) == 1:
            axis = ax
        else:
            if n == 1 or m == 1:
                axis = ax[i]
            else:
                axis = ax[i // m, i % m]
        ##start plotting
        # if len(data_list)==1:
        if to_plot == 0:
            im = axis.imshow(data_list_smoothed, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        elif to_plot < 0:
            rng = default_rng()
            if abs(to_plot) < data_list_smoothed.shape[0]:
                random_indices = rng.choice(data_list_smoothed.shape[0], size=np.abs(to_plot), replace=False)
                im = axis.imshow(data_list_smoothed[random_indices], cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
            else:
                random_indices = rng.choice(data_list_smoothed.shape[0], size=data_list_smoothed.shape[0], replace=False)
                im = axis.imshow(data_list_smoothed[random_indices], cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        elif to_plot > 0:
            im = axis.imshow(data_list_smoothed[data_list_smoothed.shape[0]//2 - to_plot//2:data_list_smoothed.shape[0]//2 + to_plot//2],
                           cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        axis.set_yticks([])
        axis.set_xticks([0, axvline, data_list_smoothed.shape[1]])
        axis.axvline(x=axvline, lw=1, ls='--', c='grey')
        axis.spines['left'].set_visible(False)
        axis.set_xticklabels([str(round(x / 60)) for x in [0, axvline - 1, (data_list_smoothed.shape[1]) - 1]])
        axis.set_title(labels[i])
        if trace_average == 'Y':
            new_ax = axis.inset_axes([0, -0.45, 1, 0.4], sharex=axis)
            mean = savgol_filter(np.nanmean(data_list_smoothed, axis=0), 61, 2)
            new_ax.plot(mean, lw=2, c='k', zorder=2.5)
            axis.spines['bottom'].set_visible(False)
            new_ax.set_xticks([0, axvline, data_list_smoothed.shape[1]])
            new_ax.set_xticklabels([str(round(x / 60)) for x in [0, axvline - 1, (data_list_smoothed.shape[1]) - 1]])
            new_ax.axvline(x=axvline, lw=1, ls='--', c='grey')
            new_ax.axhline(0, lw=1, c='grey')
            new_ax.set_ylim(top=trace_ylim[1]+20, bottom=trace_ylim[0])
            new_ax.set_yticks([trace_ylim[0], 0, trace_ylim[1]])
            new_ax.spines['bottom'].set_bounds(high=data_list_smoothed.shape[1], low=0)
            new_ax.spines['left'].set_bounds(high=trace_ylim[0], low=trace_ylim[1])
            new_ax.spines['left'].set_position(('outward', 4))
            new_ax.spines['bottom'].set_position(('outward', 2))
            if fly_indices!=None:
                each_fly_response = []
                for ids in np.unique(fly_indices[i]):
                    indices = np.argwhere(fly_indices[i] == ids).flatten()
                    if indices.shape[0]<0: ## this can be either 5 or 10, minimum number of trials for each fly
                        print('this fly {id} has too few trials'.format(id=ids))
                    else:
                        this_fly_response = savgol_filter(np.nanmean(gaussian_filter1d(data_list_smoothed[indices], sigma=3, axis=1), axis=0), 61, 2)
                        each_fly_response.append(this_fly_response)
                        new_ax.plot(this_fly_response, c='grey', lw=0.5, alpha=0.7)
                        # new_ax.plot(this_fly_response, lw=0.5, alpha=0.7)
                new_ax.plot(np.nanmean(each_fly_response, axis=0), c='grey', lw=1, alpha=0.7)
            else:
                error = sem(data_list_smoothed, axis=0)
                new_ax.fill_between(x=np.arange(0, mean.shape[0], 1), y1=mean + error, y2=mean - error, facecolor='grey', alpha=0.1)

        mean = np.nanmean(data_list_smoothed, axis=0)
        error = sem(data_list_smoothed, axis=0)
        median = np.median(data_list_smoothed, axis = 0)
        mean_trace.append(mean)
        error_trace.append(error)
        median_trace.append(median)

    i += 1
    # if n != 1 and m != 1:
    #     cbar_axis = fig.colorbar(im, cax=ax[i // m, i % m], shrink=0.1)
    # else:
    #     cbar_axis = fig.colorbar(im, cax=ax[len(data_list)], shrink=0.1)

    if n!=1 and m!=1:
        while (i < n * m):
            fig.delaxes(ax[i // m, i % m])
            i += 1
    fig.savefig(r'{}/{}'.format(Dir, filename + '.png'))
    if PDF=='Y':
        fig.savefig(r'{}/{}'.format(Dir, filename + '.pdf'), transparent=True, dpi=600)
    elif PDF == 'SVG' or PDF=='Y':
        fig.savefig(r'{}/{}'.format(Dir, filename + '.pdf'), transparent=True, dpi=600)
        fig.savefig(r'{}/{}'.format(Dir, filename + '.svg'), transparent=True, dpi=600)
    return mean_trace, error_trace, median_trace


def draw_rastor_plot_with_cluster_id(data_list, Dir, filename, cmap=ListedColormap(all_color_data['PiYG7']), vmin=0, vmax=100, labels=[],
                     smoothed='Y', normalize='N', to_plot=0, PDF='N', axvline=None, cluster_ids=[], cluster_cmap='nipy_spectral'):
    """
    Draw a rastor plot of the provided data, but, with a colormap showing flyid beside the rastor. Bad old code.
    :param data_list: The data to be plotted, presorted according to the fly ids(nested list)
    :type data_list: list of length k (k=#stimulus conditions) of list (l=#flies) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :returns: mean and median of data and error (compared to the linear prediction, relevant for angular speed data only)
    :rtype: (list of numpy.ndarray, list of numpy.ndarray, list of numpy.ndarray)
    """
    if axvline is None:
        axvline=data_list[0].shape[1]//2
    n, m = subplot_arrangement(len(data_list)) # n rows and m columns

    w_ratio = np.zeros(m*2)
    w_ratio[::2] = 1
    w_ratio[1::2] = 0.1
    fig, ax = plt.subplots(n, m*2, figsize=(figsize[0] * m, figsize[1] * n), gridspec_kw={'height_ratios': [1] * n, 'width_ratios': list(w_ratio), 'wspace':0})

    mean_trace = []
    error_trace = []
    median_trace = []
    for i in range(len(data_list)):
        if smoothed=='Y':
            data_list_smoothed = gaussian_filter1d(data_list[i], sigma=3, axis=1)
        else:
            data_list_smoothed = np.array(data_list[i])

        if normalize=='Y':
            data_list_smoothed = minmax_scale(data_list_smoothed, (-1, 1), axis = 0)
            vmax = 1
            vmin = -1
        else:
            pass

        if to_plot==0:
            pass
        elif to_plot > 0:
            data_list_smoothed = data_list_smoothed[:to_plot]
        elif to_plot < 0:
            data_list_smoothed = data_list_smoothed[to_plot:]

        l = data_list_smoothed.shape[1]
        if n == 1 or m == 1:
            im = ax[i*2].imshow(data_list_smoothed, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
            ax[i*2].set_yticks([])
            ax[i*2].set_xticks([0, data_list[i].shape[1] // 2, data_list[i].shape[1]])
            ax[i*2].set_title(labels[i])
            ax[i*2].axvline(x=axvline, lw=1, ls='--', c='grey')
            ax[(i*2)+1].imshow(np.array(cluster_ids[0][i].astype(np.float)).reshape(-1, 1), cmap=cluster_cmap, aspect='auto')
            ax[(i*2)+1].set_yticks([])
            ax[(i*2)+1].set_xticks([])
            ax[(i*2)+1].spines['left'].set_visible(False)
            ax[(i*2)+1].spines['bottom'].set_visible(False)
        else:
            im = ax[i // m, (i % m)*2].imshow(data_list_smoothed, cmap=cmap, vmin=vmin, vmax=vmax,aspect='auto')
            ax[i // m,(i % m) *2].set_yticks([])
            ax[i // m, (i % m)*2].set_xticks([0, data_list[i].shape[1] // 2, data_list[i].shape[1]])
            ax[i // m, (i % m)*2].set_title(labels[i])
            ax[i // m, (i % m)*2].axvline(x=axvline, lw=1, ls='--', c='grey')
            ax[i // m, ((i % m)*2)+1].imshow(np.array(cluster_ids[0][i].astype(np.float)).reshape(-1,1), cmap=cluster_cmap, aspect='auto')
            ax[i // m, ((i % m)*2)+1].spines['left'].set_visible(False)
            ax[i // m, ((i % m)*2)+1].spines['bottom'].set_visible(False)

        mean = np.nanmean(data_list_smoothed, axis=0)
        error = sem(data_list_smoothed, axis=0)
        median = np.median(data_list_smoothed, axis = 0)
        mean_trace.append(mean)
        error_trace.append(error)
        median_trace.append(median)

    # i += 1
    # if n != 1 and m != 1:
    #     cbar_axis = fig.colorbar(im, cax=ax[i // m, i % m], shrink=0.1)
    # else:
    #     cbar_axis = fig.colorbar(im, cax=ax[len(data_list)], shrink=0.1)
    if n!=1 and m!=1:
        while (i < n * m):
            fig.delaxes(ax[i // m, i % m])
            i += 1

    fig.savefig(r'{}/{}'.format(Dir, filename + '.png'))
    if PDF=='Y':
        fig.savefig(r'{}/{}'.format(Dir, filename + '.pdf'), transparent=True, dpi=600)
    elif PDF=='SVG':
        fig.savefig(r'{}/{}'.format(Dir, filename + '.pdf'), transparent=True, dpi=600)
        fig.savefig(r'{}/{}'.format(Dir, filename + '.svg'), transparent=True, dpi=600)

    return mean_trace, error_trace, median_trace