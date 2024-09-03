import numpy as np
from scipy.stats import sem
from scipy.stats import linregress

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import statsmodels.api as sm

import matplotlib
plt.style.use('roshan_style.mplstyle')
figsize = matplotlib.rcParams['figure.figsize']

def draw_trace_plot(data_list, Dir, filename, labels, top=100, bottom=-100, cmap='Accent', smoothed='N', highlight=None,ls=[], sem_or_not='N', normalize='N', rotate='N',
                    sort=[0, 1], draw_trace='Y', filter='Y', draw_median='Y', PDF='N', axvline=None, draw_linear_sum='N', start=0, end=-1, fly_indices=None, plot_all_together='Y',
                    draw_mean='Y', draw_var='N'):
    """
    Draw a line plots of the provided data,used for plotting speed or angular speed. Bad old code.
    Also, perform linear prediction and plot the error. 
    And return average cumulative errors, defined as the mean/sum of prediction error over some time period.
    :param data_list: The data to be plotted
    :type data_list: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :returns: fig1, ax1, mean_trace, error_trace, median_trace, cumulative_error_list, prediction_error
    :rtype: (matplotlib.figure.Figure, matplotlib.axes._subplots.AxesSubplot, list of numpy.ndarray, list of numpy.ndarray, list of numpy.ndarray, list, list)
    """
    if axvline is None:
        axvline=data_list[0].shape[1]//2
    n, m = subplot_arrangement(len(data_list))
    if not ls:
        ls = ['solid'] * len(labels)

    if cm.get_cmap(cmap).N == 256:
        colors = [cm.get_cmap(cmap)(i) for i in np.linspace(0, 1, len(data_list))]
        colors = [cm.get_cmap(cmap)(i) for i in np.linspace(0, 1, len(data_list))]
    else:
        colors = [cm.get_cmap(cmap)(i) for i in np.arange(0, len(data_list))]
    cumulative_error_list = []
    prediction_error = []
    if rotate == 'N':
        fig, ax = plt.subplots(n, m, figsize=(figsize[0]* m, figsize[1]* n + 3))
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(3, 4.5))
        mean_trace = []
        error_trace = []
        median_trace = []
        mean_each_fly_trace = []
        selected_fly_ids_for_diff_condns = []
        for i in range(len(data_list)):
            if smoothed=='Y':
                data_list_smoothed = gaussian_filter1d(data_list[i], sigma=3, axis=1)
            else:
                data_list_smoothed = np.array(data_list[i])
            if normalize=='Y':
                data_list_smoothed = minmax_scale(data_list_smoothed, (-1, 1), axis=1)
                vmax = 1
                vmin = -1
            else:
                pass

            if len(sort)==2 and sort[0]!=sort[1]:
                indices = np.argsort(np.nanmean(data_list_smoothed[:, sort[0]:sort[1]], axis=1)).flatten()
                data_list_smoothed = data_list_smoothed[indices]
                colors2 = [cm.get_cmap('summer')(i) for i in np.linspace(0, 1, data_list_smoothed.shape[0])]
            else:
                colors2 = [(0.5, 0.5, 0.5) for i in range(data_list_smoothed.shape[0])]
                pass

            if start == 0 and end == -1:
                pass
            else:
                print(data_list_smoothed.shape)
                data_list_smoothed = data_list_smoothed[:, start:end]

            l = data_list_smoothed.shape[0]
            if len(data_list)==1:
                axis = ax
            else:
                if n == 1 or m == 1:
                    axis = ax[i]
                else:
                    axis = ax[i//m, i%m]
            ## start plotting
            mean = np.nanmean(data_list_smoothed, axis=0)
            median = np.median(data_list_smoothed, axis=0)
            error = sem(data_list_smoothed, axis=0)
            if draw_trace == 'Y':
                for t, traces in enumerate(data_list_smoothed):
                    axis.plot(traces, lw=0.5, c='grey')
            else:
                pass
            percentile = int(0.2*l)
            if len(sort) == 2:
                if draw_var == 'Y':
                    axis.plot(median_abs_deviation(data_list_smoothed, axis=0), lw=2, c='k')
                else:
                    mean_top20 = np.nanmean(data_list_smoothed[: percentile], axis=0)
                    mean_bottom20 = np.nanmean(data_list_smoothed[-1 * percentile:], axis=0)
                    error_top20 = sem(data_list_smoothed[: percentile], axis=0)
                    error_bottom20 = sem(data_list_smoothed[-1 * percentile:], axis=0)
                    if filter == 'Y':
                        axis.plot(savgol_filter(mean_top20, 31, 2), lw=1.5, c='green')
                        axis.plot(savgol_filter(mean_bottom20, 31, 2), lw=1.5, c='goldenrod')
                        if sem_or_not == 'Y':
                            axis.fill_between(x=np.arange(0, mean_top20.shape[0], 1), y1=mean_top20 + error_top20, y2=mean_top20 - error_top20, facecolor='green', alpha=0.2)
                            axis.fill_between(x=np.arange(0, mean_bottom20.shape[0], 1), y1=mean_bottom20 + error_bottom20, y2=mean_bottom20 - error_bottom20, facecolor='goldenrod', alpha=0.2)
                    else:
                        axis.plot(mean_top20, lw=1.5, c='green')
                        axis.plot(mean_bottom20, lw=1.5, c='goldenrod')
                        if sem_or_not == 'Y':
                            axis.fill_between(x=np.arange(0, mean_top20.shape[0], 1), y1=mean_top20 + error, y2=mean_top20 - error, facecolor='green', alpha=0.2)
                            axis.fill_between(x=np.arange(0, mean_bottom20.shape[0], 1), y1=mean_bottom20 + error, y2=mean_bottom20 - error, facecolor='goldenrod', alpha=0.2)
            if draw_mean == 'Y':
                if filter == 'Y':
                    axis.plot(savgol_filter(mean, 31, 2), lw=3.5, c='k')
                else:
                    axis.plot(mean, lw=3.5, c='k')
                if sem_or_not=='Y':
                    axis.fill_between(x=np.arange(0, mean.shape[0], 1), y1=mean + error, y2=mean - error, facecolor='k', alpha=0.2)
            if draw_median=='Y':
                axis.plot(np.median(data_list_smoothed, axis=0), lw=2, c='k', ls='--')
            if highlight == None:
                pass
            else:
                axis.plot(np.arange(l // 2 - highlight[0], l // 2 + highlight[1]), np.nanmean(data_list_smoothed, axis=0)[l // 2 - highlight[0]:l // 2 + highlight[1]], lw=2, c='r')
            axis.set_ylim(top=top, bottom=bottom)
            axis.set(aspect=(data_list_smoothed.shape[1])/(top-bottom))
            axis.axvline(x=axvline, lw=1, ls='--', c='grey')
            axis.axhline(y=0, lw=1, ls='--', c='grey', alpha=0.5)
            axis.set_ylim(top=top+20, bottom=bottom)
            axis.set_yticks([top, 0, bottom])
            axis.set_xticks([0, axvline, data_list_smoothed.shape[1]])
            axis.set_xticklabels([0, int(axvline//60), int(data_list_smoothed.shape[1]//60)])
            axis.spines['bottom'].set_bounds(high=data_list_smoothed.shape[1], low=0)
            axis.spines['left'].set_bounds(high=top, low=bottom)
            axis.spines['left'].set_position(('outward', 4))
            axis.spines['bottom'].set_position(('outward', 2))

            if draw_linear_sum=='Y':
                pass
            else:
                if draw_var=='Y':
                    ax1.plot(median_abs_deviation(data_list_smoothed, axis=0), lw=2, c=colors[i], label=labels[i], ls=ls[i])
                else:
                    ax1.plot(mean, lw=2, c=colors[i], label=labels[i], ls=ls[i])

                ax1.plot(mean, lw=2, c=colors[i], label=labels[i], ls=ls[i])
            ax2.plot(median, lw=2, c=colors[i], label=labels[i], ls=ls[i])
            if sem_or_not == 'Y':
                ax1.fill_between(x=np.arange(0, mean.shape[0], 1), y1=mean + error, y2=mean - error, facecolor=colors[i], alpha=0.1)
                ax2.fill_between(x=np.arange(0, median.shape[0], 1), y1=median + error, y2=median - error, facecolor=colors[i], alpha=0.1)
            mean_trace.append(mean)
            error_trace.append(error)
            median_trace.append(median)
            if fly_indices != None:
                each_fly_trace = []
                selected_fly_ids = []
                for ids in np.unique(fly_indices[i]):
                    indices = np.argwhere(fly_indices[i] == ids).flatten()
                    if indices.shape[0] < 5:
                        pass
                    else:
                        each_fly_trace.append(np.nanmean(data_list_smoothed[indices], axis=0))
                        selected_fly_ids.append(ids)
                selected_fly_ids_for_diff_condns.append(selected_fly_ids)
                mean_each_fly_trace.append(each_fly_trace)
        if draw_linear_sum=='Y':
            simple_sum = mean_trace[0]+mean_trace[1]
            # ax3[0].plot(mean_trace[0], lw=1.5, c='gray', ls='dotted', label='FrontToBack', alpha=1)
            # ax3[0].plot(mean_trace[1], lw=1.5, c='gray', ls='dashed', label='BackToFront', alpha=1)
            # ax3[1].plot(simple_sum, lw=1.5, c='gray', label='FtB + BtF', alpha=0.5)
            # ax3[1].plot(mean_trace[2], lw=2, c='k', label='Full rotation')
            # ax1.plot(simple_sum, lw=1.5, c='gray', label='FtB + BtF', alpha=0.5)
            # ax1.plot(mean_trace[2], lw=2, c='k', label='Full rotation')
            # ax1.fill_between(x=np.arange(0,mean_trace[2].shape[0]), y1=simple_sum, y2=mean_trace[2], color='gainsboro', alpha=0.7)
            # ax3[1].fill_between(x=np.arange(0, mean_trace[2].shape[0]), y1=simple_sum, y2=mean_trace[2], color='gainsboro', alpha=0.7)
            prediction_error=mean_trace[2] - simple_sum

            prediction_error_list = []
            simple_sum_list = []
            FtB_list = []
            BtF_list = []
            full_rotation_list = []
            if fly_indices != None:
                final_selected_fly_ids = selected_fly_ids_for_diff_condns[0]
                for arrays in selected_fly_ids_for_diff_condns[1:]:
                    final_selected_fly_ids = np.intersect1d(final_selected_fly_ids, arrays)
                print(final_selected_fly_ids)
                for ids in final_selected_fly_ids:
                    simple_sum_each_fly = mean_each_fly_trace[0][np.argwhere(selected_fly_ids_for_diff_condns[0]==ids)[0,0]] + \
                                          mean_each_fly_trace[1][np.argwhere(selected_fly_ids_for_diff_condns[1]==ids)[0,0]]
                    difference_each_fly = mean_each_fly_trace[2][np.argwhere(selected_fly_ids_for_diff_condns[2]==ids)[0,0]] - simple_sum_each_fly
                    prediction_error_list.append(difference_each_fly)
                    simple_sum_list.append(simple_sum_each_fly)
                    full_rotation_list.append(mean_each_fly_trace[2][np.argwhere(selected_fly_ids_for_diff_condns[2]==ids)[0,0]])

                    early_response_difference = np.nanmean(difference_each_fly[stim_period - start + 10:stim_period - start + 10 + 30])
                    late_response_difference = np.nanmean(difference_each_fly[stim_period - start + 180 + 30:stim_period - start + 180 + 60])
                    total_difference = np.nanmean(difference_each_fly[stim_period-start:])

                    FtB_early_response = np.nanmean(mean_each_fly_trace[0][np.argwhere(selected_fly_ids_for_diff_condns[0]==ids)[0,0]][stim_period - start + 10:stim_period - start + 10 + 30])
                    FtB_late_response = np.nanmean(mean_each_fly_trace[0][np.argwhere(selected_fly_ids_for_diff_condns[0]==ids)[0,0]][stim_period - start + 180 + 30:stim_period - start + 180 + 60])
                    FtB_total_response = np.nanmean(mean_each_fly_trace[0][np.argwhere(selected_fly_ids_for_diff_condns[0]==ids)[0,0]][stim_period - start+10:])

                    BtF_early_response = np.nanmean(mean_each_fly_trace[1][np.argwhere(selected_fly_ids_for_diff_condns[1]==ids)[0,0]][stim_period - start + 10:stim_period - start + 10 + 30])
                    BtF_late_response = np.nanmean(mean_each_fly_trace[1][np.argwhere(selected_fly_ids_for_diff_condns[1]==ids)[0,0]][stim_period - start + 180 + 30:stim_period - start + 180 + 60])
                    BtF_total_response = np.nanmean(mean_each_fly_trace[1][np.argwhere(selected_fly_ids_for_diff_condns[1]==ids)[0,0]][stim_period - start+10:])

                    full_early_response = np.nanmean(mean_each_fly_trace[2][np.argwhere(selected_fly_ids_for_diff_condns[2]==ids)[0,0]][stim_period - start + 10:stim_period - start + 10 + 30])
                    full_late_response = np.nanmean(mean_each_fly_trace[2][np.argwhere(selected_fly_ids_for_diff_condns[2]==ids)[0,0]][stim_period - start + 180 + 30:stim_period - start + 180 + 60])
                    full_total_response = np.nanmean(mean_each_fly_trace[2][np.argwhere(selected_fly_ids_for_diff_condns[2]==ids)[0,0]][stim_period - start+10:])

                    FtB_list.append(mean_each_fly_trace[0][np.argwhere(selected_fly_ids_for_diff_condns[0]==ids)[0,0]])
                    BtF_list.append(mean_each_fly_trace[1][np.argwhere(selected_fly_ids_for_diff_condns[1]==ids)[0,0]])
                    cumulative_error_list.append([early_response_difference, late_response_difference, total_difference, FtB_early_response, FtB_late_response, FtB_total_response,
                                                  BtF_early_response, BtF_late_response, BtF_total_response, full_early_response, full_late_response, full_total_response])
                #subtract the pre-stimulus angular speed
                mean_FTB = np.subtract(np.nanmean(FtB_list, axis=0), np.nanmean(np.nanmean(FtB_list, axis=0)[:stim_period - start]))
                ax3[0].plot(mean_FTB, lw=1, c='gray', label='FrontToBack', alpha=1)
                ax3[0].fill_between(x=np.arange(FtB_list[0].shape[0]), y1=mean_FTB-sem(FtB_list, axis=0),y2=mean_FTB+sem(FtB_list, axis=0), color='gray', label='FrontToBack', alpha=0.5)
                mean_BTF = np.subtract(np.nanmean(BtF_list, axis=0), np.nanmean(np.nanmean(BtF_list, axis=0)[:stim_period - start]))
                ax3[0].plot(mean_BTF, lw=1, c='gray', label='BackToFront', alpha=1)
                ax3[0].fill_between(x=np.arange(BtF_list[0].shape[0]), y1=mean_BTF - sem(BtF_list, axis=0), y2=mean_BTF + sem(BtF_list, axis=0), color='gray',label='BackToFront', alpha=0.5)

                mean_simple_sum = np.subtract(np.nanmean(simple_sum_list, axis=0), np.nanmean(np.nanmean(simple_sum_list, axis=0)[:stim_period - start]))
                ax3[1].plot(mean_simple_sum, lw=1.5, c='gray', label='FtB + BtF', alpha=1)
                ax3[1].fill_between(x=np.arange(simple_sum_list[0].shape[0]), y1=mean_simple_sum - sem(simple_sum_list, axis=0), y2=mean_simple_sum + sem(simple_sum_list, axis=0),
                                   color='gray', label='FtB + BtF', alpha=0.5)

                mean_full_rotation = np.subtract(np.nanmean(full_rotation_list, axis=0), np.nanmean(np.nanmean(full_rotation_list, axis=0)[:stim_period - start]))
                ax3[1].plot(mean_full_rotation, lw=1.5, c='k', label='Full rotation', alpha=1)
                ax3[1].fill_between(x=np.arange(full_rotation_list[0].shape[0]), y1=mean_full_rotation - sem(full_rotation_list, axis=0), y2=mean_full_rotation + sem(full_rotation_list, axis=0),
                                    color='k', label='Full rotation', alpha=0.5)

                ax1.plot(mean_simple_sum, lw=1.5, c='gray', label='FtB + BtF', alpha=1)
                ax1.fill_between(x=np.arange(simple_sum_list[0].shape[0]), y1=mean_simple_sum - sem(simple_sum_list, axis=0), y2=mean_simple_sum + sem(simple_sum_list, axis=0),
                                    color='gray', label='FtB + BtF', alpha=0.5)
                ax1.plot(mean_full_rotation, lw=1.5, c='gray', label='Full rotation', alpha=1)
                ax1.fill_between(x=np.arange(full_rotation_list[0].shape[0]), y1=mean_full_rotation - sem(full_rotation_list, axis=0), y2=mean_full_rotation + sem(full_rotation_list, axis=0),
                                 color='gray', label='Full rotation', alpha=0.5)

                # ax3[1].plot(simple_sum, lw=1.5, c='gray', label='FtB + BtF', alpha=0.5)
                # ax3[1].plot(mean_trace[2], lw=2, c='k', label='Full rotation')
                # ax1.plot(simple_sum, lw=1.5, c='gray', label='FtB + BtF', alpha=0.5)
                # ax1.plot(mean_trace[2], lw=2, c='k', label='Full rotation')
                ax1.fill_between(x=np.arange(0, mean_trace[2].shape[0]), y1=mean_simple_sum, y2=mean_full_rotation, hatch="/////", zorder=2, alpha=0.3, facecolor='w', lw=0.5)
                ax3[1].fill_between(x=np.arange(0, mean_trace[2].shape[0]), y1=mean_simple_sum, y2=mean_full_rotation, hatch="/////", zorder=2, alpha=0.3, facecolor='w', lw=0.5)

                ax3[0].fill_betweenx(y=[bottom, top], x1=stim_period - start + 10, x2=stim_period - start + 10 + 30, color='k', alpha=0.25)
                ax3[0].fill_betweenx(y=[bottom, top], x1=stim_period - start + 180+30, x2=stim_period - start + 180 + 60, color='k', alpha=0.25)

            # slope, c, r, p, std_error = linregress(simple_sum, mean_trace[2])
            # ax1.text(110, -140, 'Simple sum r_squared :' + str(round(r**2, 3)), **dict(color='k'))
            ## linearity
            # x = np.vstack((mean_trace[0], mean_trace[1])).T
            # x = sm.add_constant(x)
            # y = mean_trace[2]
            # model = sm.OLS(y, x)
            # result = model.fit()
            # c, m1, m2 = result.params
            # ax1.plot(result.fittedvalues, lw=2, c='gray', label='Linear prediction')
            # # ax1.plot(simple_sum - mean_trace[2], lw=1, c='k', label='Difference', ls='--')
            # r_squared = result.rsquared
            # ax1.text(310, -120, 'Linearity : '+str(round(m1, 2))+' '+str(round(m2, 2))+' '+str(round(c, 2))+
            #         '\n'+'r_squared : '+str(round(r_squared, 2)), **dict(color='gray'))

        ax1.set_ylim(top=top+0.1*max(top, abs(bottom)), bottom=bottom)
        ax1.set_yticks([bottom, 0, top])
        ax1.spines['left'].set_bounds(low=bottom, high=top)
        ax1.legend(loc='upper left')
        ax1.set_xticks([0, axvline-1, (data_list_smoothed.shape[1])-1])
        ax1.spines['bottom'].set_bounds(low=0, high=(data_list_smoothed.shape[1]) - 1)
        ax1.set_xticklabels([str(round(x / 60)) for x in [0, axvline, (data_list_smoothed.shape[1]) - 1]])
        ax1.axvline(x=axvline-1, lw=1, ls='--', c='grey')
        ax1.axhline(0,0.05, 0.95, lw=1.5, c='grey', alpha=0.5)
        #
        ax2.set_ylim(top=top, bottom=bottom)
        ax2.legend(loc='upper left')
        ax2.set_xticks([0, axvline-1, (data_list_smoothed.shape[1])-1])
        ax2.spines['bottom'].set_bounds(low=0, high=(data_list_smoothed.shape[1])-1)
        ax2.set_xticklabels([str(round(x / 60)) for x in [0, axvline, (data_list_smoothed.shape[1]) - 1]])
        ax2.axvline(x=axvline-1, lw=1, ls='--', c='grey')
        # ax1.set_xticks([0, 2,(data_list_smoothed.shape[1])-1])
        # ax1.axvline(x=2, lw=1, ls='--', c='grey')
        # ax1.grid(visible=True, which='major', axis='y', color='grey', lw=0.3)
        #
        ax3[0].set_ylim(top=200, bottom=bottom-0.1*max(top, abs(bottom)))
        ax3[0].set_yticks([bottom, 0, top])
        ax3[0].spines['left'].set_bounds(low=bottom, high=top)
        ax3[0].legend(loc='upper left')
        ax3[0].set_xticks([0, axvline-1, (data_list_smoothed.shape[1])-1])
        ax3[0].spines['bottom'].set_bounds(low=0, high=(data_list_smoothed.shape[1]) - 1)
        ax3[0].set_xticklabels([str(round(x / 60)) for x in [0, axvline, (data_list_smoothed.shape[1]) - 1]])
        ax3[0].axvline(x=axvline-1, lw=1, ls='--', c='grey')
        ax3[0].axhline(0,0.05, 0.95, lw=1.5, c='grey', alpha=0.5)
        ax3[1].axvline(x=axvline-1, lw=1, ls='--', c='grey')
        ax3[1].axhline(0,0.05, 0.95, lw=1.5, c='grey', alpha=0.5)
        ax3[1].spines['left'].set_bounds(low=bottom, high=top)
        ax3[1].spines['bottom'].set_bounds(low=0, high=(data_list_smoothed.shape[1]) - 1)
        i += 1
        while (i < n * m):
            fig.delaxes(ax[i // m, i % m])
            i += 1
    elif rotate == 'Y':
        fig, ax = plt.subplots(n, m, figsize=(figsize[0]* m + figsize[0]/2, figsize[1]* n + figsize[1]/2), gridspec_kw={'height_ratios': [1] * n, 'width_ratios': [1] * m})
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        mean_trace = []
        error_trace = []
        median_trace = []
        for i in range(len(data_list)):
            if smoothed == 'Y':
                data_list_smoothed = gaussian_filter1d(data_list[i], sigma=3, axis=1)
            else:
                data_list_smoothed = np.array(data_list[i])

            if normalize == 'Y':
                data_list_smoothed = minmax_scale(data_list_smoothed, (0, 1), axis=1)
                vmax = 1
                vmin = 0
            else:
                pass

            if len(sort)==2:
                indices = np.argsort(np.nanmean(data_list_smoothed[:, sort[0]:sort[1]], axis=1)).flatten()
                data_list_smoothed = data_list_smoothed[indices]
                colors2 = [cm.get_cmap('summer')(i) for i in np.linspace(0, 1, data_list_smoothed.shape[0])]
            else:
                colors2 = [(0.5, 0.5, 0.5) for i in range(data_list_smoothed.shape[0])]
                pass

            if start == 0 and end == -1:
                pass
            else:
                data_list_smoothed = data_list_smoothed[:, start:end]

            l = data_list_smoothed.shape[0]
            if len(data_list) == 1:
                axis = ax
            else:
                if n == 1 or m == 1:
                    axis = ax[i]
                else:
                    axis = ax[i // m, i % m]
            if draw_trace == 'Y':
                for t,traces in enumerate(data_list_smoothed):
                    ax.plot(traces, np.arange(data_list_smoothed.shape[1]), lw=0.5, c='grey')
                    # ax.plot(traces, np.arange(data_list_smoothed.shape[1]), lw=0.5, c=colors2[t])
            else:
                pass
            percentile = int(0.2*l)
            if len(sort) == 2:
                if filter=='Y':
                    axis.plot(savgol_filter(np.nanmean(data_list_smoothed[: percentile], axis=0), 61, 2), np.arange(data_list_smoothed.shape[1]), lw=2, c='green')
                    axis.plot(savgol_filter(np.nanmean(data_list_smoothed[ -1 *percentile:], axis=0), 61, 2), np.arange(data_list_smoothed.shape[1]), lw=2, c='goldenrod')
                else:
                    axis.plot(np.nanmean(data_list_smoothed[: percentile], axis=0), np.arange(data_list_smoothed.shape[1]), lw=2, c='green')
                    axis.plot(np.nanmean(data_list_smoothed[ -1 *percentile:], axis=0), np.arange(data_list_smoothed.shape[1]), lw=2, c='goldenrod')
            axis.plot(np.nanmean(data_list_smoothed, axis=0), np.arange(data_list_smoothed.shape[1]), lw=2, c='k')
            if draw_median == 'Y':
                axis.plot(np.median(data_list_smoothed, axis=0), np.arange(data_list_smoothed.shape[1]), lw=2, c='k', ls='--')
            if highlight == None:
                pass
            else:
                axis.plot(np.arange(l // 2 - highlight[0], l // 2 + highlight[1]), np.nanmean(data_list_smoothed, axis=0)[l // 2 - highlight[0]:l // 2 + highlight[1]], lw=2, c='r')
            axis.set_xlim(left=top, right=bottom)

            if filter=='Y':
                mean = savgol_filter(np.nanmean(data_list_smoothed, axis=0), 61, 2)

            else:
                mean = np.nanmean(data_list_smoothed, axis=0)
            ax1.plot(mean, np.arange(data_list_smoothed.shape[1]), lw=2, c=colors[i], label=labels[i], ls=ls[i])
            median = np.nanmean(data_list_smoothed, axis=0)
            ax2.plot(median, np.arange(data_list_smoothed.shape[1]), lw=2, c=colors[i], label=labels[i], ls=ls[i])
            error = sem(data_list_smoothed, axis=0)
            if sem_or_not == 'Y':
                ax1.fill_betweenx(y=np.arange(0, mean.shape[0], 1), x1=mean + error, x2=mean - error, facecolor=colors[i], alpha=0.1)
                ax2.fill_betweenx(y=np.arange(0, median.shape[0], 1), x1=median + error, x2=median - error, facecolor=colors[i], alpha=0.1)
            mean_trace.append(mean)
            error_trace.append(error)
            median_trace.append(median)
            ax1.set_xlim(left=top, right=bottom)
            ax1.legend(loc='upper left')
            ax1.set_yticks([0, axvline, (data_list_smoothed.shape[1]) - 1])
            ax1.set_yticklabels([str(round(x / 60)) for x in [0, axvline, (data_list_smoothed.shape[1]) - 1]])
            ax1.axhline(y=axvline - 1, lw=1, ls='--', c='grey')
            ax2.set_xlim(left=top, right=bottom)
            ax2.legend(loc='upper left')
            ax2.set_yticks([0, axvline, (data_list_smoothed.shape[1]) - 1])
            ax2.set_yticklabels([str(round(x / 60)) for x in [0, axvline, (data_list_smoothed.shape[1]) - 1]])
            ax2.axhline(y=axvline - 1, lw=1, ls='--', c='grey')
            # ax1.set_xticks([0, 2,(data_list_smoothed.shape[1])-1])
            # ax1.axvline(x=2, lw=1, ls='--', c='grey')
            # ax1.grid(visible=True, which='major', axis='y', color='grey', lw=0.3)
        i += 1
        while (i < n * m):
            fig.delaxes(ax[i // m, i % m])
            i += 1

    fig.savefig(r'{}/{}'.format(Dir, filename + '.png'), dpi=300)
    if plot_all_together == 'Y':
        fig1.savefig(r'{}/{}'.format(Dir, filename + '_all.png'))
        fig2.savefig(r'{}/{}'.format(Dir, filename + '_all_median.png'))
    if rotate == 'N' and draw_linear_sum == 'Y':
        fig3.savefig(r'{}/{}'.format(Dir, filename + '_parts.png'), dpi=600)

    if PDF=='Y' or PDF=='SVG':
        fig.savefig(r'{}/{}'.format(Dir, filename + '.pdf'), transparent=True, dpi=600)
        if plot_all_together=='Y':
            fig1.savefig(r'{}/{}'.format(Dir, filename + '_all.pdf'), transparent=True, dpi=600)
            fig2.savefig(r'{}/{}'.format(Dir, filename + '_all_median.pdf'), transparent=True, dpi=600)
        if rotate == 'N' and draw_linear_sum=='Y':
            fig3.savefig(r'{}/{}'.format(Dir, filename + '_parts.pdf'), dpi=600, transparent=True)
    if PDF == 'SVG' or PDF=='Y':
        fig.savefig(r'{}/{}'.format(Dir, filename + '.svg'), transparent=True, dpi=600)
        if plot_all_together == 'Y':
            fig1.savefig(r'{}/{}'.format(Dir, filename + '_all.svg'), transparent=True, dpi=600)
            fig2.savefig(r'{}/{}'.format(Dir, filename + '_all_median.svg'), transparent=True, dpi=600)
        if rotate == 'N' and draw_linear_sum == 'Y':
            fig3.savefig(r'{}/{}'.format(Dir, filename + '_parts.svg'), dpi=600, transparent=True)

    return fig1, ax1, mean_trace, error_trace, median_trace, cumulative_error_list, prediction_error
