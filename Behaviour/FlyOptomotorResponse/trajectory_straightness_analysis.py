import pandas as pd


def straightness_analysis(saccade_peak, opto_response_speed, opto_pos, opto_response, Dir, filename,sex, p, condn_indices, ignore_values = [], variables=[]):
    ## we split the data either by the current condition for the stim periods
    indices, conditions = indices_by_condition([saccade_peak[:, x] for x in np.add(condn_indices, 1)], ignore_values=ignore_values)
    ## we split the data either by the previous condition for the non-stim periods
    indices_pre, conditions = indices_by_condition([np.insert(saccade_peak[:-1, x], 0, saccade_peak[0, x]) for x in np.add(condn_indices, 1)], ignore_values=ignore_values)
    straightness_list_stim = []
    fwd_walking_segments_list_stim = []
    ori_segments_list_stim = []
    straightness_list_pre = []
    fwd_walking_segments_list_pre = []
    mean_speed_segment_list = []
    mean_speed_segment_list_pre = []
    direction_list_stim = []
    direction_list_pre = []
    for i in range(len(indices)):
        saccades = saccade_peak[indices[i], p + 1:]
        speed = opto_response_speed[indices[i], p + 1:]
        pos_x = opto_pos[indices[i], :, 0]
        pos_y = opto_pos[indices[i], :, 1]
        ori = opto_response[indices[i]]
        fwd_walking_segments_stim, ori_segments_stim, mean_speed_segments, walking_segments_list, speed_segments_list, direction = \
            get_straight_trajectory_new(saccades[:, stim_period:], speed[:, stim_period:], pos_x[:, stim_period:], pos_y[:, stim_period:], ori[:, stim_period:], saccade_peak[indices[i], 0])
        walking_segments_list = [x[i] for x in walking_segments_list for i in range(len(x))]
        speed_segments_list = [x[i] for x in speed_segments_list for i in range(len(x))]
        fwd_walking_segments_list_stim.append(fwd_walking_segments_stim)
        ori_segments_list_stim.append(ori_segments_stim)
        direction_list_stim.append(direction)
        straightness, disp, dist = compute_straightness(fwd_walking_segments_stim)
        # straightness, total_dist, new_fwd_walking_segments = compute_straightness_new(walking_segments_list, speed_segments_list)
        straightness_list_stim.append(straightness)
        mean_speed_segment_list.append(mean_speed_segments)
        ##
        saccades = saccade_peak[indices_pre[i], p + 1:]
        speed = opto_response_speed[indices_pre[i], p + 1:]
        pos_x = opto_pos[indices_pre[i], :, 0]
        pos_y = opto_pos[indices_pre[i], :, 1]
        ori = opto_response[indices_pre[i]]
        fwd_walking_segments_pre, ori_segments, mean_speed_segments, walking_segments_list, speed_segments_list, direction = \
            draw_straight_trajectory_new(saccades[:, :stim_period], speed[:, :stim_period], pos_x[:, :stim_period], pos_y[:, :stim_period], ori[:, :stim_period], saccade_peak[indices_pre[i], 0])
        fwd_walking_segments_list_pre.append(fwd_walking_segments_pre)
        straightness, disp, dist = compute_straightness(fwd_walking_segments_pre)
        # straightness, disp, dist = compute_straightness(walking_segments_list)
        direction_list_pre.append(direction)
        straightness_list_pre.append(straightness)
        mean_speed_segment_list_pre.append(mean_speed_segments)
        # straightness, total_dist_list = compute_straightness_new(fwd_walking_segments)
        # straightness_list.append(np.average(straightness, weights=total_dist_list))
    draw_violin_plot(straightness_list_stim, make_labels(conditions), conditions, Dir, filename, ylim=[0, 20], colors=['gray']*len(straightness_list_stim))
    walking_segments_stim = draw_trajectory_straight(fwd_walking_segments_list_stim,direction_list_stim, Dir, filename='stim', labels=make_labels(conditions), color='g')
    walking_segments_pre = draw_trajectory_straight(fwd_walking_segments_list_pre,direction_list_pre, Dir, filename='pre', labels=make_labels(conditions), color='r')
    # # return fwd_walking_segments_list_stim, mean_speed_segment_list
    # fig, ax = plt.subplots(1, 1, figsize=(figsize, figsize*(3/2)), squeeze=True)
    # n = len(straightness_list_stim)
    # positions = np.linspace(0, n - 1, n)
    # result1 = ax.violinplot(straightness_list_stim, positions, showmeans=True, showextrema=False, widths=0.5)
    # labels = make_labels(conditions)
    # for i in range(len(result1['bodies'])):
    #     result1['bodies'][i].set_facecolor('g')
    # ax.set_ylim(top=15, bottom=0)
    # # ax.scatter(positions, straightness_list)
    # ax.set_xticks(np.arange(0, n))
    # ax.set_xticklabels(labels, rotation = 45, ha='right')
    # ax.set_xlabel('Contrast')
    # ax.set_ylabel('Straightness')
    # fig.tight_layout()
    # fig.savefig(r'{}/{}'.format(Dir, filename +'_straigtness.png'))
    # straightness_mean_list = [np.mean(x) for x in straightness_list_stim]
    #
    # fig, ax = plt.subplots(1, 1, figsize=(figsize, figsize*(3/2)), squeeze=True)
    # n = len(straightness_list_stim)
    # positions = np.linspace(0, n - 1, n)
    # result1 = ax.violinplot(straightness_list_pre, positions, showmeans=True, showextrema=False, widths=0.5)
    # labels = make_labels(conditions)
    # for i in range(len(result1['bodies'])):
    #     result1['bodies'][i].set_facecolor('g')
    # ax.set_ylim(top=25, bottom=0)
    # # ax.scatter(positions, straightness_list)
    # ax.set_xticks(np.arange(0, n))
    # ax.set_xticklabels(labels, rotation = 45, ha='right')
    # ax.set_xlabel('Contrast')
    # ax.set_ylabel('Straightness')
    # fig.tight_layout()
    # fig.savefig(r'{}/{}'.format(Dir, filename + sex + 'straigtness_prestim' + '.png'))
    # straightness_mean_list_prestim = [np.mean(x) for x in straightness_list_pre]

    # Dir_processed_data = os.path.join(Dir, 'ProcessedData')
    # saccade_time_data = {'variables':variables, 'conditions':conditions, 'straightness_stim':straightness_list_stim,'straightness_prestim':straightness_list_prestim}
    # np.save(os.path.join(Dir_processed_data, filename+'_straightness.npy'), saccade_time_data, allow_pickle=True)
    return straightness_list_stim, straightness_list_pre, conditions, walking_segments_stim, walking_segments_pre


def straightness_analysis_hemi(saccade_peak, opto_response_speed, opto_pos, optomotor_response, opto_response_angular_speed,orientation, Dir, filename, p, condn_indices,
                               ignore_values = [], variables=[]):
    Dir_processed_data = os.path.join(os.path.dirname(Dir), 'ProcessedData')
    if not os.path.isdir(Dir_processed_data):
        os.mkdir(Dir_processed_data)
    print(len(opto_pos))

    window = 60
    ##
    indices, rejected_indices = reject_nonresponsive_trials(opto_response_angular_speed, opto_response_speed, p=p, stim_period=stim_period)
    opto_response = optomotor_response[indices]
    opto_response_speed = opto_response_speed[indices]
    opto_response_angular_speed = opto_response_angular_speed[indices]
    saccade_peak = saccade_peak[indices]
    opto_pos = opto_pos[indices]
    orientation = orientation[indices]

    ##
    opto_response_list = np.apply_along_axis(lambda x: np.subtract(x, x[x.shape[0] // 2]), axis=1, arr=opto_response[:, p + 4:])
    ## add the expected response collumn after p collumns, defaut is 1
    opto_response, exp_type, labels = get_expected_response_hemi(opto_response, p=p, translation='N')
    # opto_response, exp_type, labels = get_expected_response_quarter(opto_response, p=p)
    # opto_response, exp_type, labels = get_expected_response_3sections(opto_response, p=p)
    ##
    saccade_peak[:, 0] = opto_response[:, p]
    indices, conditions = indices_by_condition([exp_type], ignore_values=[[5]]) ## take FtB, BtF and full_rotation experiments
    indices_pre = indices
    ##
    ## multiply by expected direction such that expected values are positive and opposite responses are negative
    opto_response_list = np.multiply(opto_response_list, opto_response[:, p].reshape(-1, 1) * -1)
    opto_response_angular_speed_list = np.multiply(opto_response_angular_speed[:, p+3:], opto_response[:, p].reshape(-1, 1) * -1)

    ## draw example trajectories
    # Dir_trajectory = os.path.join(Dir, 'all_trajectories'+labels[0])
    # if not os.path.isdir(Dir_trajectory):
    #     os.mkdir(Dir_trajectory)
    # for i in indices[0]:
    #     ## only take trials in which the fly turns atleast 1 whole rotation
    #     if opto_response_list[i, -1] < -360 and abs(opto_response_list[i, 0])<abs(opto_response_list[i, -1])/5:
    #         to_plot = i
    #         print(to_plot)
    #         draw_smooth_trajectory([opto_pos[to_plot]], [saccade_peak[to_plot][p+2:]], [opto_response_angular_speed_list[to_plot]],
    #                            [opto_response_list[to_plot]], Dir_trajectory, filename+str(to_plot), PDF='N')
    Dir_trajectory = os.path.join(Dir, 'all_trajectories'+labels[2])
    if not os.path.isdir(Dir_trajectory):
        os.mkdir(Dir_trajectory)
    print(Dir_trajectory)
    # for i in indices[2]:
    #     ## only take trials in which the fly turns atleast 1 whole rotation
    #     if opto_response_list[i, -1] > 360 and abs(opto_response_list[i, 0])<abs(opto_response_list[i, -1])/5:
    #         to_plot = i
    #         if to_plot == 597:
    #             saccade_indices = np.argwhere(saccade_peak[to_plot][p+2:] != 0).flatten()
    #             print(saccade_indices)
    #             np.set_printoptions(precision=3, suppress=True)
    #             # print(opto_response_angular_speed_list[to_plot])
    #             # saccade_peak[to_plot][p+2:][399] = opto_response_angular_speed_list[to_plot][399]
    #             saccade_peak[to_plot][p + 2:][434] = opto_response_angular_speed_list[to_plot][434]
    #             saccade_peak[to_plot][p + 2:][472] = opto_response_angular_speed_list[to_plot][472]
    #             saccade_peak[to_plot][p + 2:][524] = opto_response_angular_speed_list[to_plot][524]
    #             saccade_peak[to_plot][p + 2:][535] = opto_response_angular_speed_list[to_plot][535]
    #             saccade_peak[to_plot][p + 2:][519] = 0
    #             draw_smooth_trajectory([opto_pos[to_plot]], [saccade_peak[to_plot][p+2:]], [opto_response_angular_speed_list[to_plot]], [opto_response_list[to_plot]],
    #                                    [orientation[to_plot, p + 4:]], Dir_trajectory, filename+str(to_plot)+'arrrows', PDF='Y', arrows='Y')
    #             draw_smooth_trajectory([opto_pos[to_plot]], [saccade_peak[to_plot][p+2:]], [opto_response_angular_speed_list[to_plot]], [opto_response_list[to_plot]],
    #                                    [orientation[to_plot, p + 4:]], Dir_trajectory, filename+str(to_plot), PDF='Y', arrows='N')
    #             return opto_response_angular_speed_list[to_plot], opto_response_list[to_plot], opto_response[to_plot, p + 4:]
    # exit()
    prestim_time = 5

    straightness_list_stim = []
    fwd_walking_segments_list_stim = []
    ori_segments_list_stim = []
    straightness_list_pre = []
    fwd_walking_segments_list_pre = []
    ori_segments_list_pre = []
    mean_speed_segment_list = []
    mean_speed_segment_list_pre = []
    direction_list_stim = []
    direction_list_pre = []
    data = {}
    mean_trajectory_list_per_fly_per_stim = []
    mean_straightness_stim_list_per_fly_per_stim = []
    mean_straightness_prestim_list_per_fly_per_stim = []

    for i in range(len(indices)):
        saccades = saccade_peak[indices[i], p + 1:]
        speed = opto_response_speed[indices[i], p + 3:]
        pos_x = opto_pos[indices[i], :, 0]
        pos_y = opto_pos[indices[i], :, 1]
        print(pos_x.shape, orientation.shape)
        data[labels[i].replace(' ', '')] = {'pos_x': pos_x, 'pos_y': pos_y, 'ori': orientation[indices[i], p+3:], 'stimulus_direction': saccade_peak[indices[i], 0]}

        ori = orientation[indices[i], p + 3:]
        # fwd_walking_segments_stim, ori_segments, mean_speed_segments, walking_segments_list, speed_segments_list, direction = \
        #     get_straight_trajectory_new(saccades[:, stim_period - 5:], speed[:, stim_period - 5:], pos_x[:, stim_period - 5:], pos_y[:, stim_period - 5:], ori[:, stim_period - 5:],
        #                                 saccade_peak[indices[i], 0], window=window)
        # fwd_walking_segments_stim, ori_segments, mean_speed_segments, walking_segments_list, speed_segments_list, direction = get_straight_trajectory_new('no', speed[:, stim_period - 5:],
        #         pos_x[:, stim_period - 5:], pos_y[:, stim_period - 5:], ori[:, stim_period - 5:], saccade_peak[indices[i], 0], window=window)
        fwd_walking_segments_stim, ori_segments_stim, mean_speed_segments, walking_segments_list, speed_segments_list, direction = get_straight_trajectory_new('start', speed[:, stim_period - prestim_time:],
                pos_x[:, stim_period - prestim_time:], pos_y[:, stim_period - prestim_time:], ori[:, stim_period - prestim_time:], saccade_peak[indices[i], 0], window=window)
        ## data for each fly, comment out if not needed
        fwd_walking_segments_per_fly = []
        direction_per_fly = []
        fly_id = opto_response[:, p + 2]
        fly_id_indices = [np.array(fly_id[i]).astype('int') for i in indices]
        mean_trajectory_list_per_fly = []
        mean_straightness_per_fly = []
        mean_straightness_prestim_per_fly = []
        for ids in np.unique(fly_id):
            fly_indices = np.argwhere(fly_id_indices[i] == ids).flatten()
            if fly_indices.shape[0] == 0:
                mean_straightness_per_fly.append(np.nan)
                mean_straightness_prestim_per_fly.append(np.nan)
                continue
            else:
                pass
            try:
                fwd_walking_segments_stim_per_fly, _, _, _, _, _ = get_straight_trajectory_new(saccades[:, stim_period - prestim_time:], speed[fly_indices, stim_period - prestim_time:], pos_x[fly_indices, stim_period - prestim_time:],
                                                    pos_y[fly_indices, stim_period - prestim_time:], ori[fly_indices, stim_period - prestim_time:], saccade_peak[indices[i], 0][fly_indices], window=window)
            except:
                fwd_walking_segments_stim_per_fly = np.empty((0, 1))
            if fwd_walking_segments_stim_per_fly.shape[0]!=0:
                straightness_stim_per_fly, _, _ = compute_straightness(fwd_walking_segments_stim_per_fly)
                mean_straightness_per_fly.append(np.mean(straightness_stim_per_fly))
            else:
                mean_straightness_per_fly.append(np.nan)
            try:
                fwd_walking_segments_prestim_per_fly, _, _, _, _, _ = get_straight_trajectory_new(saccades[:, :stim_period - prestim_time], speed[fly_indices, :stim_period - prestim_time], pos_x[fly_indices, :stim_period - prestim_time],
                                                    pos_y[fly_indices, :stim_period - prestim_time], ori[fly_indices, :stim_period - prestim_time], saccade_peak[indices[i], 0][fly_indices], window=window)
            except:
                fwd_walking_segments_prestim_per_fly = np.empty((0, 1))
            if fwd_walking_segments_prestim_per_fly.shape[0]!=0:
                straightness_prestim_per_fly, _, _ = compute_straightness(fwd_walking_segments_prestim_per_fly)
                mean_straightness_prestim_per_fly.append(np.mean(straightness_prestim_per_fly))
            else:
                mean_straightness_prestim_per_fly.append(np.nan)

            fwd_walking_segments_stim_per_fly, ori_segments_per_fly, mean_speed_segments_per_fly, walking_segments_list_per_fly, speed_segments_lis_per_flyt, direction_per_fly = \
                get_straight_trajectory_new('start', speed[fly_indices, stim_period - prestim_time:], pos_x[fly_indices, stim_period - prestim_time:], pos_y[fly_indices, stim_period - prestim_time:], ori[fly_indices, stim_period - prestim_time:],
                                            saccade_peak[indices[i], 0][fly_indices], window=window)
            if len(fwd_walking_segments_stim_per_fly)!=0:
                trajectory_per_fly = draw_trajectory_straight([fwd_walking_segments_stim_per_fly], [direction_per_fly], [ori_segments_per_fly], Dir, filename='trajectory_per_fly_', labels=labels, color='g', save='N')
                mean_trajectory_list_per_fly.append([np.mean(trajectory_per_fly[0][:, 0, :], axis=0), np.mean(trajectory_per_fly[0][:, 1, :], axis=0)])
        mean_trajectory_list_per_fly = np.array(mean_trajectory_list_per_fly)
        mean_trajectory_list_per_fly_per_stim.append(mean_trajectory_list_per_fly)
        #
        mean_straightness_per_fly=np.array(mean_straightness_per_fly)
        mean_straightness_prestim_per_fly=np.array(mean_straightness_prestim_per_fly)
        mean_straightness_stim_list_per_fly_per_stim.append(mean_straightness_per_fly)
        mean_straightness_prestim_list_per_fly_per_stim.append(mean_straightness_prestim_per_fly)
        ##
        walking_segments_list = [x[i] for x in walking_segments_list for i in range(len(x))]
        speed_segments_list = [x[i] for x in speed_segments_list for i in range(len(x))]
        fwd_walking_segments_list_stim.append(fwd_walking_segments_stim)
        ori_segments_list_stim.append(ori_segments_stim)
        direction_list_stim.append(direction)
        straightness, disp, dist = compute_straightness(fwd_walking_segments_stim)
        # straightness, total_dist, new_fwd_walking_segments = compute_straightness_new(walking_segments_list, speed_segments_list)
        straightness_list_stim.append(straightness)
        mean_speed_segment_list.append(mean_speed_segments)
        ##
        saccades = saccade_peak[indices_pre[i], p + 1:]
        speed = opto_response_speed[indices_pre[i], p + 3:]
        pos_x = opto_pos[indices_pre[i], :, 0]
        pos_y = opto_pos[indices_pre[i], :, 1]
        ori = orientation[indices_pre[i]]
        # fwd_walking_segments_pre, ori_segments, mean_speed_segments, walking_segments_list, speed_segments_list, direction = \
        #     get_straight_trajectory_new(saccades[:, :stim_period - 5], speed[:, :stim_period - 5], pos_x[:, :stim_period - 5], pos_y[:, :stim_period - 5], ori[:, :stim_period - 5], saccade_peak[indices_pre[i], 0], window=window)
        # fwd_walking_segments_pre, ori_segments, mean_speed_segments, walking_segments_list, speed_segments_list, direction = \
        #     get_straight_trajectory_new('no', speed[:, :stim_period - 5], pos_x[:, :stim_period - 5], pos_y[:, :stim_period - 5], ori[:, :stim_period - 5], saccade_peak[indices_pre[i], 0], window=window)
        fwd_walking_segments_pre, ori_segments_pre, mean_speed_segments, walking_segments_list, speed_segments_list, direction = \
            get_straight_trajectory_new('start', speed[:, :stim_period - prestim_time], pos_x[:, :stim_period - prestim_time], pos_y[:, :stim_period - prestim_time], ori[:, :stim_period - prestim_time],
                                        saccade_peak[indices_pre[i], 0], window=window)
        fwd_walking_segments_list_pre.append(fwd_walking_segments_pre)
        ori_segments_list_pre.append(ori_segments_pre)
        # straightness, disp, dist = compute_straightness(walking_segments_list)
        direction_list_pre.append(direction)
        straightness_list_pre.append(straightness)
        mean_speed_segment_list_pre.append(mean_speed_segments)
        # straightness, total_dist_list = compute_straightness_new(fwd_walking_segments)
        # straightness_list.append(np.average(straightness, weights=total_dist_list))
    scipy.io.savemat(os.path.join(Dir, filename + 'fly_positions.mat'), data)
    print("saved")
    return 1
    mean_straightness_prestim_list_per_fly_per_stim = np.nanmean(mean_straightness_prestim_list_per_fly_per_stim, axis=0)
    mean_straightness_stim_list_per_fly_per_stim = np.array(mean_straightness_stim_list_per_fly_per_stim)

    ##
    # genotype = filename.split('_hemi_')[0]
    # if genotype=='CantonS_18':
    #     genotype = 'CantonS'
    # print(genotype)
    # df = pd.DataFrame({'straightness': mean_straightness_prestim_list_per_fly_per_stim, 'stim_type': ['pre']*mean_straightness_prestim_list_per_fly_per_stim.shape[0],
    #                    'genotype': [genotype]*mean_straightness_prestim_list_per_fly_per_stim.shape[0]})
    # labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    # for i in range(mean_straightness_stim_list_per_fly_per_stim.shape[0]):
    #     df1 = pd.DataFrame({'straightness': mean_straightness_stim_list_per_fly_per_stim[i, :], 'stim_type': [labels[i]] * mean_straightness_stim_list_per_fly_per_stim.shape[1],
    #                         'genotype': [genotype]*mean_straightness_prestim_list_per_fly_per_stim.shape[0]})
    #     df = pd.concat([df, df1], ignore_index=True)
    # df.to_csv(os.path.join(Dir_processed_data, filename + '_straightness.csv'))


    mean_trajectory_list_per_fly_per_stim = np.array(mean_trajectory_list_per_fly_per_stim)
    n, m = subplot_arrangement(len(indices))
    fig, ax = plt.subplots(n, m, figsize=(m * 3, n * 3))
    for i in range(mean_trajectory_list_per_fly_per_stim.shape[0]):
        for j in range(len(mean_trajectory_list_per_fly_per_stim[i])):
            ax[i].plot(mean_trajectory_list_per_fly_per_stim[i][j, 0, :], mean_trajectory_list_per_fly_per_stim[i][j, 1, :], c='grey', lw=1, alpha=0.3)
        ax[i].plot(np.mean(mean_trajectory_list_per_fly_per_stim[i][:, 0, :], axis=0), np.mean(mean_trajectory_list_per_fly_per_stim[i][:, 1, :], axis=0), c='k', lw=2, alpha=0.8)
        ax[i].set_xlim(left=-100, right=100)
        ax[i].set_ylim(top=150, bottom=0)
        ax[i].set_xticks([])
        ax[i].set_yticks([])
        ax[i].set_title(labels[i])
        ax[i].set_aspect('equal')
    fig.savefig(os.path.join(Dir, filename+'_trajectory_per_fly.png'))
    fig.savefig(os.path.join(Dir, filename+'_trajectory_per_fly.pdf'), dpi=600, transparent=True)
    # draw_violin_plot([straightness_list_stim], labels, conditions[0][0], Dir, filename, ylim=[0, 20], colors=['gray']*len(straightness_list_stim), cluster='Y', PDF='N')
    walking_segments_stim = draw_trajectory_straight(fwd_walking_segments_list_stim, direction_list_stim, ori_segments_list_stim,
                                                     Dir, filename=filename+'stim_start', labels=labels, color='g')
    walking_segments_pre = draw_trajectory_straight(fwd_walking_segments_list_pre, direction_list_pre, ori_segments_list_pre,
                                                    Dir, filename=filename+'pre_start', labels=labels, color='r')
    # return 1
    fig, ax = plt.subplots(1, 1, figsize=(figsize[0], figsize[1]), squeeze=True)
    n = len(straightness_list_stim)
    positions = np.linspace(0, n - 1, n)
    result1 = ax.violinplot(straightness_list_stim, positions, showmeans=True, showextrema=False, widths=0.5)
    for i in range(len(result1['bodies'])):
        result1['bodies'][i].set_facecolor('g')
        result1['bodies'][i].set_alpha(0.5)
    result2 = ax.violinplot(straightness_list_pre, positions, showmeans=True, showextrema=False, widths=0.5)
    for i in range(len(result2['bodies'])):
        result2['bodies'][i].set_facecolor('gray')
        result2['bodies'][i].set_alpha(0.5)
    ax.set_ylim(top=1, bottom=0.5)
    # ax.scatter(positions, straightness_list)
    ax.set_xticks(np.arange(0, n))
    ax.set_xticklabels(labels)
    # ax.set_xlabel('Contrast')
    ax.set_ylabel('Straightness')
    fig.tight_layout()
    fig.savefig(r'{}/{}'.format(Dir, filename +'_straigtness.png'))
    straightness_mean_list = [np.mean(x) for x in straightness_list_stim]

    # Dir_processed_data = os.path.join(Dir, 'ProcessedData')
    # saccade_time_data = {'variables':variables, 'conditions':conditions, 'straightness_stim':straightness_list_stim,'straightness_prestim':straightness_list_prestim}
    # np.save(os.path.join(Dir_processed_data, filename+'_straightness.npy'), saccade_time_data, allow_pickle=True)
    return straightness_list_stim, straightness_list_pre, conditions, walking_segments_stim, walking_segments_pre


def get_straight_trajectory_new(saccades, speed, pos_x, pos_y, ori, direction, window=30):
    if type(saccades)!=str:
        """generates forward walking segments of size 20 frames for further analysis"""
        indices = np.argwhere(saccades)
        walking_segments_list = []
        speed_segments_list = []
        ori_segments_list = []
        """collect all the walking segments around saccades"""
        for k in range(speed.shape[0]):
            """each k is one trial"""
            j = 0
            walking_segments = []
            speed_segments = []
            ori_segments = []
            for i in indices[np.argwhere(indices[:, 0] == k).flatten(), 1]:
                if i > j + 25:
                    """i = position of the saccades
                    j = beginning of window, the window ends 5 frames before the beginning of the saccade"""
                    walking_segments.append([pos_x[k, j:i - 5], pos_y[k, j:i - 5]])
                    speed_segments.append([speed[k, j:i - 5]])
                    ori_segments.append([ori[k, j:i - 5]])
                j = i + 10
            walking_segments_list.append(walking_segments)
            speed_segments_list.append(speed_segments)
            ori_segments_list.append(ori_segments)
    else:
        if saccades=='no':
            walking_segments_list = []
            speed_segments_list = []
            ori_segments_list = []
            for k in range(speed.shape[0]):
                """each k is one trial"""
                walking_segments = [[pos_x[k, :], pos_y[k, :]]]
                speed_segments = [[speed[k, :]]]
                ori_segments = [[ori[k, :]]]
                walking_segments_list.append(walking_segments)
                speed_segments_list.append(speed_segments)
                ori_segments_list.append(ori_segments)
        elif saccades=='start':
            # the trajectory starts when the stimulus starts
            start=0
            walking_segments_list = []
            speed_segments_list = []
            ori_segments_list = []
            """collect all the walking segments around saccades"""
            for k in range(speed.shape[0]):
                """each k is one trial"""
                walking_segments = [[pos_x[k, start:start+window+1], pos_y[k, start:start+window+1]]]
                speed_segments = [[speed[k, start:start+window+1]]]
                ori_segments = [[ori[k, start:start+window+1]]]
                walking_segments_list.append(walking_segments)
                speed_segments_list.append(speed_segments)
                ori_segments_list.append(ori_segments)
    ##walking_segments_list is a list of size speed.shape[0] with lists of varying sizes,
    ##each with a list of size 2 (x and y) comprising of np arrays
    # walking_segments_new_list = []
    # speed_segments_new_list = []
    # ori_segments_new_list = []
    walking_segments_total_list = []
    speed_segments_total_list = []
    ori_segments_total_list = []
    direction_list = []
    for k in range(speed.shape[0]):
        for i in range(len(walking_segments_list[k])):
            for j in range(walking_segments_list[k][i][0].shape[0] // window): ## break down the segments into smaller segments of size "window"
                walking_segments_total_list.append([walking_segments_list[k][i][0][j * window:(j + 1) * window], walking_segments_list[k][i][1][j * window:(j + 1) * window]])
                speed_segments_total_list.append(speed_segments_list[k][i][0][j * window:(j + 1) * window])
                ori_segments_total_list.append(ori_segments_list[k][i][0][j * window:(j + 1) * window])
                direction_list.append(direction[k])
    speed_segments_total_list = np.array(speed_segments_total_list)
    walking_segments_total_list = np.array(walking_segments_total_list)
    ori_segments_total_list = np.array(ori_segments_total_list)
    direction_list = np.array(direction_list)
    mean_speed_segments = np.mean(speed_segments_total_list, axis=1)
    distance = np.apply_along_axis(lambda x: np.diff(x), axis=2, arr=walking_segments_total_list)

    # index of fwd_segments where the mean speed more than 5 and there are no jumps
    indices = np.setdiff1d(np.argwhere(mean_speed_segments > 2).flatten(), np.unique(np.argwhere(np.abs(distance) > 20)[:, 0])) ## choose only segment where the fly is running and NOT jumping
    fwd_walking_segments = walking_segments_total_list[indices]
    ori_segments = ori_segments_total_list[indices]
    direction_list = direction_list[indices]
    return fwd_walking_segments, ori_segments, mean_speed_segments[indices], walking_segments_list, speed_segments_list, direction_list


def compute_straightness_new(fwd_walking_segments):
    fwd_segments=[]
    for i in range(len(fwd_walking_segments)):
        fwd_segments += fwd_walking_segments[i]
    straightness = []
    total_dist_list = []
    for i in range(len(fwd_segments)):
        per_dist = []
        dist = []
        dist_x = np.diff(fwd_segments[i][0])
        dist_y = np.diff(fwd_segments[i][1])
        total_dist = np.sum(np.sqrt(np.square(dist_x) + np.square(dist_y)))
        for j in range(10, fwd_segments[i][0].shape[0]-10):
            p1 = np.array([fwd_segments[i][0][j-10], fwd_segments[i][1][j-10]])
            p2 = np.array([fwd_segments[i][0][j+10], fwd_segments[i][1][j+10]])
            p3 = np.array([fwd_segments[i][0][j], fwd_segments[i][1][j]])
            per_dist.append(abs(np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1)))
            dist.append(abs(np.linalg.norm(p1-p2)))
        if sum(per_dist)==0:
            pass
        else:
            straightness.append(sum(dist)/sum(per_dist))
            total_dist_list.append(total_dist)
    return straightness, total_dist_list


def compute_straightness(fwd_segments):
    """compute the how straight a given forward segment is by dividing the length of the segment by the perpendicular distance at the midpoint;
    the larger the value, the straighter the segment"""
    dist_x = np.apply_along_axis(lambda x: np.diff(x), axis=1, arr=fwd_segments[:, 0, :])
    dist_y = np.apply_along_axis(lambda x: np.diff(x), axis=1, arr=fwd_segments[:, 1, :])
    dist = np.sum(np.sqrt(np.square(dist_x) + np.square(dist_y)), axis=1)
    disp_x = np.apply_along_axis(lambda x: np.subtract(x[0], x[-1]), axis=1, arr=fwd_segments[:, 0, :])
    disp_y = np.apply_along_axis(lambda x: np.subtract(x[0], x[-1]), axis=1, arr=fwd_segments[:, 1, :])
    disp = np.sqrt(np.square(disp_x) + np.square(disp_y))
    straightness = np.divide(np.abs(disp), np.abs(dist)) ##maximum value is 1
    # straightness = np.divide(np.multiply(np.subtract(np.abs(dist), np.abs(disp)), np.abs(disp)), np.abs(dist)) ##maximum value is 1
    # non_nan_indices = [~np.isnan(straightness)]
    # print(non_nan_indices)
    # print()
    straightness = straightness[~np.isnan(straightness)]
    disp = disp[~np.isnan(straightness)]
    dist = dist[~np.isnan(straightness)]
    return straightness, disp, dist
    # p1 = np.array(np.stack((fwd_segments[:,0,0], fwd_segments[:,1,0]), axis = 1))
    # p2 = np.array(np.stack((fwd_segments[:,0,-1], fwd_segments[:,1,-1]), axis = 1))
    # p3 = np.array(np.stack((fwd_segments[:,0, fwd_segments.shape[2]//2], fwd_segments[:,1,fwd_segments.shape[2]//2]), axis =1))
    # per_dist = abs(np.linalg.norm(np.cross(p2-p1, p1-p3))/np.linalg.norm(p2-p1))
    # straightness = np.divide(np.abs(dist), np.abs(per_dist))
    # return straightness, per_dist, dist
