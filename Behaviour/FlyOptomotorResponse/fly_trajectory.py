import os
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib import cm

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def draw_fly_trajectory(pos_x, pos_y, Dir, filename='', radius1=510, radius2=510):
    """
    Create and save a plot of the fly trajectory
    :param pos_x: x value of fly position
    :type pos_x: numpy.ndarray
    :param pos_y: y value of fly position
    :type pos_y: numpy.ndarray
    """
    fig, ax = plt.subplots()
    ax.plot(pos_x, pos_y, lw=0.4, color='k')
    circle = mpatches.Circle((512, 512), radius=radius1, color='gray', fill=False, lw=2)
    circle2 = mpatches.Circle((512, 512), radius=radius2, color='k', fill=False, lw=2)
    ax.add_artist(circle2)
    ax.add_artist(circle)
    ax.set_frame_on(False)
    ax.set_xlim(left=0, right=1050)
    ax.set_ylim(bottom=0, top=1050)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(os.path.join(Dir, filename + '_trajectory.png'))
    return 1


def draw_trajectory_straight(fwd_walking_segments_list, direction_list, ori_segments_list, Dir, labels, filename='', color='r', save='Y'):
    """
    Create a plot of trajectories that have been straightened.
    :param fwd_walking_segments_list: The fly trajectories.
    :type fwd_walking_segments_list: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, 2, m) where n is the number of trials and m is the number of frames.
    :param direction_list: The direction of the stimulus.
    :type direction_list: ist of length k (k=#stimulus conditions) of numpy.ndarray of shape (n,) where n is the number of trials
    :param ori_segments_list: The orientation of the fly corresponding to the trajectories.
    :type ori_segments_list: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :returns: straigtened trajectories (such that they point in the same direction at the start)
    :rtype: list of length k (k=#stimulus conditions) of numpy.ndarray of shape (n, 2, m) where n is the number of trials and m is the number of frames.
    """
    def cart2pol(x, y):
        theta = np.arctan2(y, x)
        rho = np.hypot(x, y)
        return theta, rho

    def pol2cart(theta, rho, offset=None):
        if offset is None:
            offset = [0, 0]
        x = rho * np.cos((theta / 180) * np.pi)
        y = rho * np.sin((theta / 180) * np.pi)
        return x - offset[0], y - offset[1]
    n, m = subplot_arrangement(len(fwd_walking_segments_list))
    if n==1:
        n=2
    if m==1:
        m=2
    fig, ax = plt.subplots(n, m, figsize=(m*6, n*6))
    transformed_walking_segments = []
    for j in range(len(fwd_walking_segments_list)):
        # print('drawing trajectories...' + str(fwd_walking_segments_list[j].shape))
        fwd_walking_segments = np.apply_along_axis(lambda x: np.subtract(x, x[0]), axis=2, arr=fwd_walking_segments_list[j])
        ori_segments = ori_segments_list[j]
        pos_x_list = []
        pos_y_list = []
        theta_list = []
        rho_list = []
        if len(direction_list)!=0:
            direction = direction_list[j]
            for i in range(fwd_walking_segments.shape[0]):
                theta_bout, rho_bout = cart2pol(fwd_walking_segments[i, 0, :], fwd_walking_segments[i, 1, :])
                theta_bout = np.multiply((theta_bout / np.pi), 180)
                offset_angle, _ = cart2pol(fwd_walking_segments[i, 0, 5], fwd_walking_segments[i, 1, 5])
                offset_angle = ((offset_angle / np.pi) * 180) - 90
                # offset_angle = np.median(ori_segments[i, :5])
                theta_bout = np.subtract(theta_bout, offset_angle) % 360
                theta_bout = ((theta_bout - 180) % 360) - 180
                pos_x_relative, pos_y_relative = pol2cart(theta_bout, rho_bout)
                # flip only the x-axis depending on the direction of stimulus, this is a bit weird to do
                # we are pretending that for cw stimulus, the animal should move to the right
                if direction[i] == -1:
                    pass
                else:
                    # print('1')
                    pos_x_relative = pos_x_relative * -1
                final_theta_bout, final_rho_bout = cart2pol(pos_x_relative, pos_y_relative)
                final_theta_bout = np.multiply((final_theta_bout / np.pi), 180)
                final_theta_bout = ((final_theta_bout - 180) % 360) - 180
                if n==1 and m==1:
                    if np.abs(np.amax(pos_x_relative)) > 150 or np.abs(np.amax(pos_y_relative)) > 150:
                        pass
                    else:
                        ax.plot(pos_x_relative, pos_y_relative, c='grey', lw=0.3, alpha=0.5)
                else:
                    if np.abs(np.amax(pos_x_relative)) > 200 or np.abs(np.amax(pos_y_relative)) > 200:
                        pass
                    else:
                        ax[j // m, j % m].plot(pos_x_relative, pos_y_relative, c='grey', lw=0.3, alpha=0.5)
                pos_x_list.append(pos_x_relative)
                pos_y_list.append(pos_y_relative)
                # theta_bout = theta_bout%360
                # print(np.max(theta_bout), np.min(theta_bout))
                theta_list.append(final_theta_bout)
                rho_list.append(final_rho_bout)
                # print('here')
        else:
            for i in range(fwd_walking_segments.shape[0]):
                theta_bout, rho_bout = cart2pol(fwd_walking_segments[i, 0, :], fwd_walking_segments[i, 1, :])
                theta_bout = np.multiply((theta_bout / np.pi), 180)
                offset_angle, _ = cart2pol(fwd_walking_segments[i, 0, 5], fwd_walking_segments[i, 1, 5])
                offset_angle = ((offset_angle / np.pi) * 180) - 90
                # offset_angle = np.median(ori_segments[i, :5])
                theta_bout = np.subtract(theta_bout, offset_angle) % 360
                pos_x_relative, pos_y_relative = pol2cart(theta_bout, rho_bout)
                if np.abs(np.amax(pos_x_relative))>150 or np.abs(np.amax(pos_y_relative))>150:
                    pass
                else:
                    ax[j//m, j % m].plot(pos_x_relative, pos_y_relative, c=color, lw=0.3, alpha=0.5)
                pos_x_list.append(pos_x_relative)
                pos_y_list.append(pos_y_relative)
                theta_list.append(theta_bout)
                rho_list.append(rho_bout)
        if n ==1 and m==1:
            # mean_pos_x, mean_pos_y = pol2cart(np.nanmean(theta_list, axis=0), np.nanmean(rho_list, axis=0))
            ax.plot(np.nanmean(pos_x_list, axis=0), np.nanmean(pos_y_list, axis=0), c='k', lw=2, solid_capstyle='round')
            # ax.plot(mean_pos_x, mean_pos_y, c='g', lw=2, solid_capstyle='round')
            ax.set_xlim(left=-200, right=200)
            ax.set_ylim(top =200, bottom =-50)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(labels[j])
            ax.set_aspect('equal')
        else:
            # mean_pos_x, mean_pos_y = pol2cart(np.nanmean(theta_list, axis=0), np.nanmean(rho_list, axis=0))
            ##
            diff_x = np.sum(np.diff(pos_x_list, axis=1), axis=0)
            diff_y = np.sum(np.diff(pos_y_list, axis=1), axis=0)
            new_pos_x = np.cumsum(diff_x)
            new_pos_y = np.cumsum(diff_y)
            ax[j // m, j % m].plot(new_pos_x, new_pos_y, c='r', lw=2, solid_capstyle='round')
            ##
            ax[j // m, j % m].plot(np.nanmean(pos_x_list, axis=0), np.nanmean(pos_y_list, axis=0), c='k', lw=2, solid_capstyle='round')
            # ax[j // m, j % m].plot(mean_pos_x, mean_pos_y, c='g', lw=2, solid_capstyle='round')
            ax[j // m, j % m].set_xlim(left=-200, right=200)
            ax[j // m, j % m].set_ylim(top=200, bottom=-50)
            # ax[j // m, j % m].set_xticks([])
            # ax[j // m, j % m].set_yticks([])
            ax[j // m, j % m].set_title(labels[j])
            ax[j // m, j % m].set_aspect('equal')
        transformed_walking_segments.append(np.stack((pos_x_list, pos_y_list), axis=1))
    plt.subplots_adjust(left=1, right=1)
    j += 1
    while (j < n * m):
        fig.delaxes(ax[j // m, j % m])
        j += 1
    if save=='Y':
        fig.savefig(os.path.join(Dir, filename+'fwd_segments.png'))
        fig.savefig(os.path.join(Dir, filename +'fwd_segments.pdf'), dpi=600, transparent=True)
    else:
        pass
    return transformed_walking_segments


def draw_smooth_trajectory(pos_data, saccade_peak, angular_speed, ori, ori_original, Dir, filename='', colors=[], PDF='Y', arrows='N'):
    """
    draw example fly paths, with angular speed and total change in orientation, len(data) gives how many paths are drawn
    :param pos_data: list of [[pos_x],[pos_y]]
    :param Dir:
    :param filename:
    :param colors:
    :return:
    """
    def cart2pol(x, y):
        theta = np.arctan2(y, x)
        rho = np.hypot(x, y)
        return theta, rho

    def pol2cart(theta, rho, offset=None):
        if offset is None:
            offset = [0, 0]
        x = rho * np.cos((theta / 180) * np.pi)
        y = rho * np.sin((theta / 180) * np.pi)
        return x - offset[0], y - offset[1]
    fig, ax = plt.subplots(1,2, figsize=(8, 3), gridspec_kw={'width_ratios':[1,2]})
    ax_ori = ax[1].twinx() ## duplicate the axis for angular speed
    ax[1].set_ylabel('Angular speed (degrees/s)', labelpad=1)
    ax_ori.set_ylabel('Cumulative turns (degrees)', horizontalalignment='right', fontdict={'color':'grey'}, labelpad=1)
    start = 0
    stim=300
    end = 600
    cmap=ListedColormap(all_color_data['RdGy8']) ## starts with red ends with blue
    colors=[cm.get_cmap(cmap)(i) for i in np.linspace(1,0, end-start)]
    for i in range(len(pos_data)):
        val = np.linspace(0, 1, end-start)
        fwd_walking_segments = np.subtract(pos_data[i][start:end], pos_data[i][start]).T
        theta_bout, rho_bout = cart2pol(fwd_walking_segments[0, :], fwd_walking_segments[1, :])
        theta_bout = np.multiply((theta_bout / np.pi), 180)
        # theta_bout = np.subtract(theta_bout, ori_segments_list[i, 0] % 360)
        # theta_bout = np.subtract(theta_bout, ori_segments_list[i, 0] % 360)
        offset_angle, _ = cart2pol(fwd_walking_segments[0, 10], fwd_walking_segments[1, 10])
        offset_angle = ((offset_angle / np.pi) * 180) - 90
        theta_bout = np.subtract(theta_bout, offset_angle)
        pos_x_relative, pos_y_relative = pol2cart(theta_bout, rho_bout)

        saccade_indices = np.argwhere(saccade_peak[i] != 0).flatten()
        for j in range(1, pos_x_relative.shape[0]):
            if arrows=='Y':
                ax[0].plot(pos_x_relative[j - 1:j + 1], pos_y_relative[j - 1:j + 1], color=colors[j], lw=0.5, solid_capstyle='round', zorder=2, alpha=0.2)
            else:
                ax[0].plot(pos_x_relative[j - 1:j + 1], pos_y_relative[j - 1:j + 1], color=colors[j], lw=1.5, solid_capstyle='round', zorder=2, alpha=1)
            # ax[1].plot([j - 1, j], angular_speed[i][j - 1:j + 1], color=colors[j], lw=1.5, solid_capstyle='round')
            # ax[1].plot([j - 1, j], angular_speed[i][j - 1:j + 1], color=colors[j], lw=1.5, solid_capstyle='round')
        # to mark the beginning of the stimulus
        ax[0].scatter(pos_x_relative[stim], pos_y_relative[stim], marker='x', s=100, c='k', zorder=5)
        ax[1].plot(angular_speed[i][start:end], color='k', lw=1.5)
        ax_ori.plot(ori[i][start:end], color='grey', lw=2.5, alpha=0.5)

        cmap = make_color_pastel(ListedColormap(all_color_data['PiYG5']))

        for s_i in saccade_indices:
            if angular_speed[i][s_i]>0:
                color=cmap(0.95)
            else:
                color=cmap(0.05)
            ax[1].scatter(s_i, angular_speed[i][s_i], color=color, lw=1, marker='o', alpha=0.8, linewidths=0)
            if arrows=='Y':
                ax[0].plot(pos_x_relative[s_i-8:s_i], pos_y_relative[s_i-8:s_i], color=color, lw=0.5, solid_capstyle='round', zorder=2, alpha=0.2)
                # ax[0].scatter(pos_x_relative[s_i-5], pos_y_relative[s_i-5], color='slateblue', marker='o', alpha=0.3, zorder=2, s=10)
            else:
                ax[0].plot(pos_x_relative[s_i-8:s_i], pos_y_relative[s_i-8:s_i], color=color, lw=1.5, solid_capstyle='round', zorder=2, alpha=0.75)

        com = [np.nanmean(pos_x_relative), np.nanmean(pos_y_relative)] # center of mass, subtract to keep the trajectory more or less centered
        max_x = np.amax(np.subtract(pos_x_relative, com[0]))+5
        max_y = np.amax(np.subtract(pos_y_relative, com[1]))+5
        min_x = np.amin(np.subtract(pos_x_relative, com[0]))-5
        min_y = np.amin(np.subtract(pos_y_relative, com[1]))-5
        if arrows == 'Y':
            ax[0].quiver(pos_x_relative[::2], pos_y_relative[::2], np.cos(((ori_original[i][::2]-offset_angle+90)/180)*np.pi), np.sin(((ori_original[i][::2]-offset_angle+90)/180)*np.pi),
                         color='grey', angles='xy', scale_units='xy', scale=0.05, pivot='mid', alpha=0.7)
        else:
            pass
    ax[0].set_xlim(left=com[0]+min_x, right=com[0]+max_x)
    # ax[0].set_ylim(top=com[1]+max_y, bottom=com[1]+min_y)
    ax[0].set_ylim(top=com[1]+max_y, bottom=com[1]+min_y)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_frame_on(False)
    ax[0].set(aspect=1)
    ax[0].plot([com[0]+max_x-100, com[0]+max_x-10], [com[1]+min_y,com[1]+min_y], lw=3, c='k')

    max_response = np.amax(np.abs(angular_speed))
    max_response = ((max_response//50)+1)*50
    ax[1].set_ylim(top=max_response+50, bottom=-max_response-50)
    # ax[1].set_frame_on(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['left'].set_bounds(low=-max_response, high=max_response)
    ax[1].set_yticks([max_response, -max_response])
    ax[1].set_xticks([])
    ax[1].axvline(x=300, lw=2, c='grey', ls='--', alpha=0.3)
    ax[1].axhline(y=0, lw=2, c='grey', alpha=0.3)

    max_response = np.amax(np.abs(ori))
    max_response = ((max_response//50)+1)*50
    ax_ori.set_ylim(top=max_response+10, bottom=-max_response-10)
    ax_ori.set_xticks([])
    ax_ori.set_yticks([max_response, -max_response])
    ax_ori.set_frame_on(False)
    # sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('PuOr'), norm = matplotlib.colors.Normalize(-1, 1))
    # plt.colorbar(sm)
    fig.savefig(os.path.join(Dir, filename + '_scatter_trajectorey.png'), dpi=600)
    if PDF=='Y':
        fig.savefig(os.path.join(Dir, filename + '_scatter_trajectorey.pdf'), transparent=True, dpi=600)
    return fig, ax