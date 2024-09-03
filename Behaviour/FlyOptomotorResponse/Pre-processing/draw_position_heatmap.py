import numpy as np
from scipy.ndimage import gaussian_filter

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def draw_heatmap_2d(data, color_map = 'viridis', vmin = 0, vmax=0, Dir = r'D:', file = 'figure.png', radius = 512):
    data = gaussian_filter(data, sigma=1)
    fig, ax = plt.subplots(1, 1, squeeze=True)
    im = plt.imshow(data, origin='lower', extent=[0, 2*radius, 0, 2*radius], cmap=color_map, vmin = vmin, vmax=vmax)
    circle1 = mpatches.Ellipse((radius, radius), width=10, height=40, color='silver', fill=True)
    circle2 = mpatches.Ellipse((radius, radius+15), width=10, height=10, color='silver', fill=True)
    clip = mpatches.Ellipse((radius, radius), width=2*radius, height=2*radius, transform=ax.transData)
    circle3 = mpatches.Ellipse((radius, radius), width=1.5* radius, height=1.5*radius, color='k', fill=False)
    ax.add_artist(circle1)
    ax.add_artist(circle2)
    im.set_clip_path(clip)
    for pos in ['left','right','top','bottom']:
        ax.spines[pos].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    # cax = fig.add_axes(ax, sharey=True)
    # cax.set_xticks([vmin, vmax])
    # cax.set_yticklabels(["{:.2e}".format(vmin), "{:.2e}".format(vmax)])
    cb = fig.colorbar(im, ticks=[vmin, vmax], format='%.2e', aspect=15, pad=0.05)
    cb.set_ticks([])
    fig.savefig(r'{}/{}'.format(Dir, file))
    return fig


def draw_position_heatmap(list_x, list_y, bins_rad = 50, bins_theta = 50, radius = 512, vmin=None, vmax=None, color_map = 'viridis', Dir = r'D:/', file= 'figure.png'):
    heatmap, xedges, yedges = np.histogram2d(list_x, list_y, bins=(bins_theta, bins_rad), range = [[0,2*radius],[0,2*radius]], density=True)
    # np.sum(heatmap)
    heatmap_new = heatmap
    fig = draw_heatmap_2d(heatmap, color_map=color_map, vmin=vmin, vmax=vmax, Dir=Dir, file=file, radius = radius)
    return fig, heatmap_new