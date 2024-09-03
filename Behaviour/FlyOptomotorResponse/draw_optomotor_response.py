def draw_optomotor_response(opto_response, angular_speed, speed, Dir, indices=None, filename='', labels=[], top=500, bottom=-500, cmap='Accent', ls=[], rotate='N', stim_period=600, PDF='N'):
    """
    Create and save many plots related to analysis of optomotor response
    :param opto_response: The orientatio data
    :type opto_response: numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param angular_speed: The angular speed data.
    :type angular_speed: numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param speed: The speed data.
    :type speed: numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param Dir: The directory to save the plots to.
    :type Dir: str
    :param indices: The indices of the trials that belong to each stimulus condition. 
    :type indices: list of length k(where k is #stimulus conditions) of 1-d numpy.ndarray of varying lengths
    :param stim_period: The period after which the rotation starts, in frames, 5s at 60Hz is 300 frames.
    :rotate: 'Y' if the data needs to be rotated by 90 degree(i.e. time is along y-axis), 'N' otherwise.
    :type rotate: str
    :returns: Figure, axis, mean orientation response, standard error of mean orientation response, 
    median orientation response, mean angular speed, standard error of mean angular speed, median angular speed
    :rtype: (matplotlib.figure.Figure, matplotlib.axis.Axis, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray)
    """
    if indices==None:
        indices = [0]
        opto_response_list = [opto_response]
        angular_speed_list = [angular_speed]
        speed_list = [speed]
    else:
        opto_response_list = []
        angular_speed_list = []
        speed_list = []
        for index in indices:
            opto_response_list.append(opto_response[index])
            angular_speed_list.append(angular_speed[index])
            speed_list.append(speed[index])
    if not ls:
        ls = ['solid'] * len(labels)
    fig, ax = plt.subplots()
    for i in range(len(indices)):
        angular_speed_list_smoothed = gaussian_filter1d(angular_speed, sigma=3, axis=1)
        im = ax.imshow(angular_speed_list_smoothed, cmap=ListedColormap(all_color_data['RdBu3']), vmin=-250, vmax=250, aspect='auto')
    cbar_axis = fig.colorbar(im)
    ax.set_yticks([])
    fig.savefig(r'{}/{}'.format(Dir, filename+'_total_angular_speed'+'.png'))

    # draw_rastor_plot(angular_speed_list, Dir, filename + '_angular_speed_list_unsorted', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
    #                  vmin=-200, vmax=200, labels=labels, smoothed='Y', normalize='N', sort=[], trace_average='Y', start=180, end=-1, axvline=120, PDF=PDF)
    fig2, ax2, opto_response_mean, opto_response_sem, opto_response_median, _, _ = draw_trace_plot(opto_response_list, Dir, filename+'_opto_response', labels=labels, top=1000, bottom=-1000, draw_trace='N',
                                        cmap=cmap, smoothed='N', highlight=None, ls=ls, sem_or_not='Y', normalize='N', rotate=rotate, sort=[stim_period-30, stim_period], draw_median = 'N', start=180, end=-1, axvline=120, PDF=PDF)
    _, _, angular_speed_mean, angular_speed_sem, angular_speed_median, _, _ = draw_trace_plot(angular_speed_list, Dir, filename + '_angular_speed', labels=labels, top=top, bottom=bottom, cmap=cmap,
                            smoothed='Y', highlight=None, ls=ls, sem_or_not='N', normalize='N', rotate='N', sort=[stim_period-30, stim_period], draw_trace='N', draw_median = 'N', start=180, end=-1, axvline=120, PDF=PDF)
    draw_trace_plot(angular_speed_list, Dir, filename + '_angular_speed_unsorted', labels=labels, top=top, bottom=bottom, cmap=cmap, smoothed='Y', highlight=None, ls=ls, sem_or_not='Y',
                    normalize='N', rotate='N', sort=[0], draw_trace='Y', draw_median = 'N', PDF=PDF, start=180, end=-1, axvline=120, filter='Y')
    # draw_trace_plot([x[:, stim_period-120:stim_period+120] for x in angular_speed_list], Dir, filename + '_angular_speed_immediate', labels=labels, top=top, bottom=bottom, cmap=cmap,
    #                 smoothed='Y', highlight=None, ls=ls, sem_or_not='Y', normalize='N', rotate='N', sort=[0], draw_median = 'N', filter='N', PDF=PDF)
    ## sort by speed before the trial
    opto_response_list_sorted = []
    angular_speed_list_sorted = []
    for i in range(len(opto_response_list)):
        sorted_indices = np.argsort(np.nanmean(speed_list[i][:,stim_period-30: stim_period], axis = 1)).flatten()
        opto_response_list_sorted.append(opto_response_list[i][sorted_indices])
        angular_speed_list_sorted.append(angular_speed_list[i][sorted_indices])
    draw_rastor_plot(angular_speed_list_sorted, Dir, filename + '_angular_speed_list_speed_sorted_rastor', cmap=make_color_pastel(ListedColormap(all_color_data['PiYG5']), c=0.35),
                     vmin=-200, vmax=200, labels=labels, smoothed='Y', normalize='N', sort=[], PDF=PDF, trace_average='Y', start=180, end=-1, axvline=120)
    draw_trace_plot(angular_speed_list_sorted, Dir, filename + '_angular_speed_list_speed_sorted', labels=labels, top=top, bottom=bottom, draw_median='N',
                    cmap=cmap, smoothed='Y', highlight=None, ls=ls, sem_or_not='Y', normalize='N', rotate=rotate, sort=[0,0], draw_trace='N', start=180, end=-1, axvline=120, filter='Y', PDF=PDF)
    draw_trace_plot(opto_response_list_sorted, Dir, filename + '_opto_response_speed_sorted',labels=labels, top=1000, bottom=-1000,
                    cmap=cmap, smoothed='Y', highlight=None, ls=ls, sem_or_not='Y', normalize='N', rotate=rotate, sort=[0,0], start=180, end=-1, axvline=120, PDF=PDF)
    return fig2, ax2, opto_response_mean, opto_response_sem, opto_response_median, angular_speed_mean, angular_speed_sem, angular_speed_median
