import numpy as np
import os

import scipy.signal
import matplotlib.pyplot as plt


def auto_correlation_velocities(Dir, angular_speed, side_velocity, fwd_velocity, labels, filename):
    color_list = ['salmon', 'violet', 'yellowgreen']
    fig, ax = plt.subplots(1, 3, figsize=(6, 2))
    corr_angular_speed_prestim_list = []
    corr_side_vel_prestim_list = []
    corr_fwd_vel_prestim_list = []
    for i in range(len(angular_speed)):
        corr_angular_speed_list = []
        corr_side_vel_list = []
        corr_fwd_vel_list = []
        for j in range(angular_speed[i].shape[0]):
            corr_angular_speed = scipy.signal.correlate(angular_speed[i][j, 300:], angular_speed[i][j, 300:], mode='same')
            corr_angular_speed_list.append(corr_angular_speed/np.max(np.abs(corr_angular_speed)))
            corr_side_vel = scipy.signal.correlate(side_velocity[i][j, 300:], side_velocity[i][j, 300:], mode='same')
            corr_side_vel_list.append(corr_side_vel/np.max(np.abs(corr_side_vel)))
            corr_fwd_vel = scipy.signal.correlate(fwd_velocity[i][j, 300:], fwd_velocity[i][j, 300:], mode='same')
            corr_fwd_vel_list.append(corr_fwd_vel/np.max(np.abs(corr_fwd_vel)))

            corr_angular_speed = scipy.signal.correlate(angular_speed[i][j, :300], angular_speed[i][j, :300], mode='same')
            corr_angular_speed_prestim_list.append(corr_angular_speed/np.max(np.abs(corr_angular_speed)))
            corr_side_vel = scipy.signal.correlate(side_velocity[i][j, :300], side_velocity[i][j, :300], mode='same')
            corr_side_vel_prestim_list.append(corr_side_vel/np.max(np.abs(corr_side_vel)))
            corr_fwd_vel = scipy.signal.correlate(fwd_velocity[i][j, :300], fwd_velocity[i][j, :300], mode='same')
            corr_fwd_vel_prestim_list.append(corr_fwd_vel/np.max(np.abs(corr_fwd_vel)))

        ax[0].plot(np.nanmean(corr_angular_speed_list, axis=0), label=labels[i], c=color_list[i], lw=1.5, alpha=0.6)
        ax[1].plot(np.nanmean(corr_side_vel_list, axis=0), label=labels[i], c=color_list[i], lw=1.5, alpha=0.6)
        ax[2].plot(np.nanmean(corr_fwd_vel_list, axis=0), label=labels[i], c=color_list[i], lw=1.5, alpha=0.6)
    ax[0].plot(np.nanmean(corr_angular_speed_prestim_list, axis=0), label='Prestim', c='grey', lw=1.5, alpha=0.6)
    ax[1].plot(np.nanmean(corr_side_vel_prestim_list, axis=0), label='Prestim', c='grey', lw=1.5, alpha=0.6)
    ax[2].plot(np.nanmean(corr_fwd_vel_prestim_list, axis=0), label='Prestim', c='grey', lw=1.5, alpha=0.6)
    ax[1].legend()
    # ax[0].set_ylim([-0.1, 0.2])
    # ax[1].set_ylim([-0.1, 1])
    fig.savefig(os.path.join(Dir, filename + '_auto_correlation.png'))
    return 1


def cross_correlation_between_velocities(Dir, angular_speed, side_velocity, fwd_velocity, labels, filename):
    import scipy.signal
    color_list = ['salmon', 'violet', 'yellowgreen']
    fig, ax = plt.subplots(1, 2, figsize=(4, 2))
    corr_angular_fwd_vel_prestim = []
    corr_angular_side_vel_prestim = []
    for i in range(len(angular_speed)):
        corr_angular_side_vel = []
        corr_angular_fwd_vel = []
        for j in range(angular_speed[i].shape[0]):
            corr_side = scipy.signal.correlate(angular_speed[i][j, 300:], side_velocity[i][j, 300:], mode='same')
            corr_angular_side_vel.append(corr_side/np.max(np.abs(corr_side)))
            corr_fwd = scipy.signal.correlate(np.abs(angular_speed[i][j, 300:]), fwd_velocity[i][j, 300:], mode='same')
            corr_angular_fwd_vel.append(corr_fwd/np.max(np.abs(corr_fwd)))
            corr_side = scipy.signal.correlate(angular_speed[i][j, :300], side_velocity[i][j, :300], mode='same')
            corr_angular_side_vel_prestim.append(corr_side/np.max(np.abs(corr_side)))
            corr_fwd = scipy.signal.correlate(np.abs(angular_speed[i][j, :300]), fwd_velocity[i][j, :300], mode='same')
            corr_angular_fwd_vel_prestim.append(corr_fwd/np.max(np.abs(corr_fwd)))
            # spearmanr, _ = scipy.stats.spearmanr(angular_speed[i], side_velocity[i])
            # corr_angular_side_vel.append(spearmanr)
        ax[0].plot(np.nanmean(corr_angular_side_vel, axis=0), label=labels[i], c=color_list[i], lw=1.5, alpha=0.6)
        ax[1].plot(np.nanmean(corr_angular_fwd_vel, axis=0), label=labels[i], c=color_list[i], lw=1.5, alpha=0.6)
    ax[0].plot(np.nanmean(corr_angular_side_vel_prestim, axis=0), label='Prestim', c='grey', lw=1.5, alpha=0.6)
    ax[1].plot(np.nanmean(corr_angular_fwd_vel_prestim, axis=0), label='Prestim', c='grey', lw=1.5, alpha=0.6)
    ax[1].legend()
    ax[0].set_ylim([-0.1, 0.2])
    ax[1].set_ylim([-0.1, 1])
    fig.savefig(os.path.join(Dir, filename + '_cross_correlation.png'))
    return 1
