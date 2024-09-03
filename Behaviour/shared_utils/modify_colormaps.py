import numpy as np
from matplotlib.colors import ListedColormap

import matplotlib.pyplot as plt
import matplotlib.cm as cm

def make_new_colormap_alternate(n, cmap = ['Reds', 'Blues'], repeats = 1, min=0.1):
    """
    Create a new colormap by combining two colormaps such that the colors are repeated alternately.
    """
    colors = np.zeros((n*len(cmap),4))
    for j in range(len(cmap)):
        if cm.get_cmap(cmap[j]).N == 256:
            colors_new = [cm.get_cmap(cmap[j])(i) for i in np.linspace(min, 1, n)]
        else:
            colors_new = [cm.get_cmap(cmap[j])(i) for i in np.arange(0, n)]
        colors[j::2] = colors_new
    colors = colors * repeats
    return ListedColormap(colors)


def make_new_colormap(n, cmap = ['Reds', 'Blues'], repeats = 1, min=0.25):
    """
    Create a new colormap by combining two colormaps such that the colors one after another.
    """
    colors = []
    for j in range(len(cmap)):
        if cm.get_cmap(cmap[j]).N == 256:
            colors_new = [cm.get_cmap(cmap[j])(i) for i in np.linspace(min, 1, n)]
        else:
            colors_new = [cm.get_cmap(cmap[j])(i) for i in np.arange(0, n)]
        colors = colors + colors_new
    colors = colors * repeats
    return ListedColormap(colors)


def make_color_pastel(colormap, c=0.25, n=256):
    """
    Create a pastel version of a colormap by reducing the saturation, i.e., add more white to the colors.
    :param colormap: The colormap to be pastelized.
    :type colormap: matplotlib.colors.ListedColormap
    :param c: The amount of white to be added to the colors.
    :type c: float
    : returns: The pastelized colormap.
    :rtype: matplotlib.colors.ListedColormap
    """
    colormap = (1. - c) * colormap(np.linspace(0., 1., n)) + c * np.ones((n, 4))
    return ListedColormap(colormap)


def create_colormap_and_linestyle(data, condn_indices=[2,1,0], ignore_values = []):
    """
    Create a new colormap and linestyles for plotting the data based on the conditions.
    Not an elegan way of doing this.
    :param data: The data to be plotted.
    :type data: numpy.ndarray
    :param condn_indices: The indices of the conditions in the data.
    :type condn_indices: list
    :returns: The new colormap and linestyles.
    :rtype: (matplotlib.colors.ListedColormap, list)
    """
    if not ignore_values:
        ignore_values = [[None]] * len(condn_indices)
    # colors = ['YlGn', 'PuBu','Greens', 'Reds', 'Greys', 'Blues']
    colors = ['Greens', 'Reds', 'Greys', 'Blues', 'YlGn', 'PuBu']
    linestyles = ['solid', 'dashed', 'dotted']
    if len(condn_indices) == 1:
        unique_values = list(np.unique(data[:, condn_indices[0]]))
        unique_values = [x for x in unique_values if x not in ignore_values[0]]
        # unique_values = [x for x in unique_values if x !=ignore_values[0]]
        l = len(unique_values)
        newcmp = make_new_colormap(l, [colors[0]], repeats=1)
        ls = ['solid'] * l
        return newcmp, ls
    elif len(condn_indices) == 2:
        unique_values = list(np.unique(data[:, condn_indices[0]]))
        unique_values = [x for x in unique_values if x not in ignore_values[0]]
        # unique_values = [x for x in unique_values if x != ignore_values[0]]
        l = len(unique_values)
        unique_values = list(np.unique(data[:, condn_indices[1]]))
        unique_values = [x for x in unique_values if x not in ignore_values[1]]
        # unique_values = [x for x in unique_values if x != ignore_values[1]]
        l2 = len(unique_values)
        newcmp = make_new_colormap(l2, colors[:l], repeats=1) # use l colors and divide each color into l2 portions
        ls = ['solid'] * (l*l2)
        return newcmp, ls
    elif len(condn_indices) == 3:
        unique_values = list(np.unique(data[:, condn_indices[0]]))
        unique_values = [x for x in unique_values if x not in ignore_values[0]]
        # unique_values = [x for x in unique_values if x != ignore_values[0]]
        l = len(unique_values)
        unique_values = list(np.unique(data[:, condn_indices[1]]))
        unique_values = [x for x in unique_values if x not in ignore_values[1]]
        # unique_values = [x for x in unique_values if x != ignore_values[1]]
        l2 = len(unique_values)
        unique_values = list(np.unique(data[:, condn_indices[2]]))
        unique_values = [x for x in unique_values if x not in ignore_values[2]]
        # unique_values = [x for x in unique_values if x != ignore_values[2]]
        l3 = len(unique_values)
        newcmp = make_new_colormap(l, colors[:l2], repeats=l3)
        ls = []
        for i in range(l):
            ls += [linestyles[i]] * (l3 * l2)
        return newcmp, ls
    elif len(condn_indices) > 3:
        unique_values = list(np.unique(data[:, condn_indices[0]]))
        unique_values = [x for x in unique_values if x not in ignore_values[0]]
        # unique_values = [x for x in unique_values if x != ignore_values[0]]
        l = len(unique_values)
        unique_values = list(np.unique(data[:, condn_indices[1]]))
        unique_values = [x for x in unique_values if x not in ignore_values[1]]
        # unique_values = [x for x in unique_values if x != ignore_values[1]]
        l2 = len(unique_values)
        unique_values = list(np.unique(data[:, condn_indices[2]]))
        unique_values = [x for x in unique_values if x not in ignore_values[2]]
        # unique_values = [x for x in unique_values if x != ignore_values[2]]
        l3 = len(unique_values)
        unique_values = list(np.unique(data[:, condn_indices[3]]))
        unique_values = [x for x in unique_values if x not in ignore_values[3]]
        # unique_values = [x for x in unique_values if x != ignore_values[2]]
        l4 = len(unique_values)
        newcmp = make_new_colormap(l, colors[:l2], repeats=l3*l4)
        ls = []
        for i in range(l):
            ls += [linestyles[i]] * (l3 * l2*l4)
        return newcmp, ls