import random
import json
import matplotlib
matplotlib.use('agg')
from numpy.random import default_rng
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from sklearn.preprocessing import minmax_scale
from skimage import filters
import matplotlib.patches as mpatches

import numpy as np
import math
# from colour import choose_random_color
import copy
import pandas as pd
import scipy as sps
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.collections import LineCollection
import cv2
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
# import skvideo.io as sk
import os
import math
from collections import Counter
from scipy.signal import butter, lfilter
from scipy.signal import find_peaks_cwt
from scipy.interpolate import UnivariateSpline
from scipy.stats import sem
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import median_abs_deviation
from scipy.signal import savgol_filter
import statsmodels.api as sm
import seaborn as sns
import glob
import Optomotor_data_extractor
plt.style.use('roshan_style.mplstyle')
figsize = matplotlib.rcParams['figure.figsize']
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
from mat4py import loadmat
from scipy import integrate

global stim_period
stim_period = 300

global color_dict
color_dict = json.load(open('color_dict.json'))

global all_color_data
all_color_data = loadmat('colorList.mat')['all_data']
# all_color_data = {}


def fly_close_to_wall(r, limit=512, threshold=30, period=90):
    """
    Checks if the fly is close to the wall and returns frames where it is.

    :param r: The distance of the fly from the wall.
    :type r: float
    :param limit: Radius of the arena.
    :type limit: float
    :param threshold: The distance from the wall that is considered close.
    :type threshold: float
    :param period: The number of frames to consider.
    :type period: int
    :returns: Array of frames where the fly is close to the wall.
    :rtype: numpy.ndarray
    """
    indices = np.zeros(r.shape[0]//period)
    for i in range(r.shape[0]//period):
        if np.any(r[i * period:(i + 1) * period] > (limit - threshold)):
            indices[i] = 1
    return indices


def draw_scatter_plot(list_x, list_y, Dir, filename='',labels=[], cmap='Accent', ylim=[100, -100], xlim=[50, 0], regression='Y'):
    """
    Create and save a scatter plot of the data.
    :param list_x: The x-axis data.
    :type list_x: list of numpy.ndarray
    :param list_y: The y-axis data.
    :type list_y: list of numpy.ndarray
    :regression: 'Y' if the regression line needs to be plotted, 'N' otherwise.
    :type regression: str
    """
    n, m = subplot_arrangement(len(list_x))
    fig, ax = plt.subplots(n, m, figsize=(figsize[0]* m, figsize[1]* n + 3), gridspec_kw={'height_ratios': [1] * n, 'width_ratios': [1] * m})
    fig1, ax1 = plt.subplots()
    if cm.get_cmap(cmap).N == 256:
        colors = [cm.get_cmap(cmap)(i) for i in np.linspace(0, 1, len(list_x))]
    else:
        colors = [cm.get_cmap(cmap)(i) for i in np.arange(0, len(list_y))]
    R2=0
    if len(list_x)==1:
        for i in range(len(list_x)):
            if regression == 'Y':
                R2, R2_shuffled, coeff, intercept=regression_analysis([np.array([list_x[i]]).T],[np.array([list_y[i]]).T])
                y = np.dot(np.array([[xlim[0]], [xlim[1]]]),coeff[0]) + intercept[0]
                y1 = y[0, 0]
                y2 = y[1, 0]
                ax.plot([xlim[0], xlim[1]], [y1, y2], lw=2, color='grey')
                ax1.plot([xlim[0], xlim[1]], [y1, y2], lw=2, color=colors[i], label=labels)
                ax.text(xlim[1]-((xlim[1]-xlim[0])*0.1), ylim[1]-((ylim[1]-ylim[0])*0.1), 'R2 ' + "%.3f" % R2[0])
            ax.scatter(list_x[i], list_y[i], s=20, color='grey', alpha=0.4)
            ax.set_ylim(top=ylim[0], bottom=ylim[1])
            ax.set_xlim(right=xlim[0], left=xlim[1])
            ax.set_title(labels[i])
            ax1.scatter(list_x[i], list_y[i], s=50, c=colors[i], alpha=0.4, label=labels[i]+'  R2'+"%.3f" % R2[0])
    elif n==1 or m==1:
        for i in range(len(list_x)):
            if regression == 'Y':
                R2, R2_shuffled, coeff, intercept=regression_analysis([np.array([list_x[i]]).T],[np.array([list_y[i]]).T])
                y = np.dot(np.array([[xlim[0]], [xlim[1]]]),coeff[0]) + intercept[0]
                y1 = y[0, 0]
                y2 = y[1, 0]
                ax[i].plot([xlim[0], xlim[1]], [y1, y2], lw=2, color='grey')
                ax1.plot([xlim[0], xlim[1]], [y1, y2], lw=2, color=colors[i], label=labels)
                ax[i].text(xlim[1]-((xlim[1]-xlim[0])*0.1), ylim[1]-((ylim[1]-ylim[0])*0.1), 'R2 ' +"%.3f" % R2[0])
            ax[i].scatter(list_x[i], list_y[i], s=20, c='grey', alpha=0.4)
            ax[i].set_ylim(top=ylim[0], bottom=ylim[1])
            ax[i].set_xlim(right=xlim[0], left=xlim[1])
            ax[i].set_title(labels[i])
            ax1.scatter(list_x[i], list_y[i], s=50, c=colors[i], alpha=0.4, label=labels[i]+'  R2'+"%.3f" % R2[0])
    else:
        for i in range(len(list_x)):
            if regression == 'Y':
                R2, R2_shuffled, coeff, intercept=regression_analysis([np.array([list_x[i]]).T],[np.array([list_y[i]]).T])
                y = np.dot(np.array([[xlim[0]], [xlim[1]]]),coeff[0]) + intercept[0]
                y1 = y[0, 0]
                y2 = y[1, 0]
                ax[i // m, i % m].plot([xlim[0], xlim[1]], [y1, y2], lw=2, color='grey')
                ax1.plot([xlim[0], xlim[1]], [y1, y2], lw=2, color=colors[i], label=labels)
                ax[i // m, i % m].text(xlim[1]-((xlim[1]-xlim[0])*0.1), ylim[1]-((ylim[1]-ylim[0])*0.1), 'R2 ' +"%.3f" % R2[0])
            ax[i // m, i % m].scatter(list_x[i], list_y[i], s=20, c='grey', alpha=0.4)
            ax[i // m, i % m].set_ylim(top=ylim[0], bottom=ylim[1])
            ax[i // m, i % m].set_xlim(right=xlim[0], left=xlim[1])
            ax[i // m, i % m].set_title(labels[i])
            ax1.scatter(list_x[i], list_y[i], s=50, c=colors[i], alpha=0.4, label=labels[i]+'  R2'+"%.3f" % R2[0])
    ax1.set_ylim(top=ylim[0], bottom=ylim[1])
    ax1.set_xlim(right=xlim[0], left=xlim[1])
    # ax1.legend()
    fig.savefig(os.path.join(Dir, filename + '.png'))
    fig1.savefig(os.path.join(Dir, filename + 'together.png'))
    return 1


def sort_list_of_lists(my_list, rule_list1, rule_list2):
    my_list_new = copy.deepcopy(list(my_list))
    for i in range(len(my_list_new)):
        #my_list_new[i]=np.insert(my_list_new[i],0,np.nanmean(my_list_new[i][0:60]))
        my_list_new[i] = np.insert(my_list_new[i], 0, rule_list1[i])
        my_list_new[i] = np.insert(my_list_new[i], 0, rule_list2[i])
        my_list_new[i] = list(my_list_new[i])

    def sorting_rule(val):
        return val[:2]

    my_list_new.sort(key=sorting_rule)
    my_list_new = np.array(my_list_new)
    return my_list_new


def plot_inst_stim_resp(data_list): ## pre-treated data
    list_mean = []
    n = len(data_list)
    fig,ax = plt.subplots(1,n)
    for k in range(n):
        mean_list = []
        for i in range(len(data_list[k])):
            ax[k].plot(data_list[k][i][2:], color='grey', lw=0.25)
            mean_list.append(data_list[k][i][2:])
        mean = np.nanmean(mean_list, axis=0)
        list_mean.append(mean)
        ax[k].plot(mean, color='k', lw=1)
        ax[k].set_ylim(-300,400)
        #plt.axvline(60, ls='--', lw=2)
    plt.show()
    return 1


def rotateImage(image, angle):
    """
    Rotate an image
    :param image: input image
    :type image: np.array
    :param angle: angle of rotation
    :type angle: float
    """
    image_center = tuple(np.array(image.shape[1::-1]) / 2)
    rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
    result = cv2.warpAffine(image, rot_mat, image.shape[1::-1])
    return result


def turn_bias_wrt_wall(ori, r, theta, turn):
    ## computes if the turns made by the fly are biased wrt the wall and center
    if turn == 0:
        return 0
    ori = (ori/180) * np.pi
    body_vector = [10*np.cos(ori), 10*np.sin(ori)]
    center_vector = [r*np.cos(theta), r*np.sin(theta)]
    cross = np.cross(body_vector,center_vector)
    angle = np.arcsin(float(cross)/(r*10))
    angle = (angle/np.pi) * 180

    return turn*((turn/angle)/abs(turn/angle))


def pol2cart(r,theta):
    theta = (theta/180)*np.pi
    x = abs(r)*np.cos(theta)
    y = abs(r)*np.sin(theta)
    return x, y


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(x, y)
    return(rho, phi)


def circular_transformation(x_list, y_list, ori): ## has to be corrected ori
    ori = (ori/180) * np.pi
    r = []
    for i in range(len(x_list)):
        r.append(np.sqrt(x_list[i]**2 + y_list[i]**2))
    new_x = []
    new_y = []
    dir = np.arctan2(x_list, y_list)
    for i in range(len(r)):
        new_x.append(r[i] * np.sin(dir[i]-ori))
        new_y.append(r[i] * np.cos(dir[i]-ori))

    return np.array([np.array(new_x), np.array(new_y)])


def escape_response(frame_change_looms, pos_x, pos_y, ori): ## has to be corrected ori
    escape_trace = []
    ori_trace = []
    for i in frame_change_looms:
        if i +4 +60 > len(pos_x):
            escape_trace.append(np.array([np.zeros(60),np.zeros(60)]))
        else:
            escape_trace.append(circular_transformation(np.subtract(pos_x[i+4:i+4+60], pos_x[i+4]), np.subtract(pos_y[i+4:i+4+60], pos_y[i+4]), ori[i+4]))
            ori_trace.append(np.subtract(ori[i+4:i+4+60], ori[i+4]))

    return escape_trace, ori_trace


def fix_ori_angles(ori):
    ori_new = []
    fly_ori_last = 0
    fly_frame_ori_last = 0
    for i in range(len(ori)):
        fly_ori = ori[i]
        fly_turn = fly_ori - fly_ori_last

        if fly_turn > 90:
            fly_turn = -(180 - fly_turn)
        elif fly_turn < -90:
            fly_turn = 180 + fly_turn

        fly_frame_ori = fly_frame_ori_last + fly_turn
        ori_new.append(fly_frame_ori)
        fly_ori_last = fly_ori
        fly_frame_ori_last = fly_frame_ori
    return ori_new


def regression_analysis(data_X, data_Y):
    from sklearn.linear_model import Ridge
    from sklearn.linear_model import LinearRegression as LR
    from sklearn.model_selection import train_test_split

    R2_list = []
    R2_shuffled_list = []
    coefficient_list = []
    intercept_list = []
    for i in range(len(data_X)):
        non_nan_data = ~(np.isnan(data_X[i][:, 0]) | np.isnan(data_Y[i][:, 0]) | np.isinf(data_X[i][:, 0]) | np.isinf(data_Y[i][:, 0]))
        data_X[i] = data_X[i][non_nan_data]
        data_Y[i] = data_Y[i][non_nan_data]

        X_train, X_test, y_train, y_test = train_test_split(data_X[i], data_Y[i], test_size=0.20)
        y_train_shuffled = copy.deepcopy(y_train)
        np.random.shuffle(y_train_shuffled)
        model = LR().fit(X_train, y_train)
        predicted_values = model.predict(X_test)
        R2 = model.score(X_test, y_test)

        model2 = LR().fit(X_train, y_train_shuffled)
        predicted_values2 = model2.predict(X_test)
        R2_shuffled = model2.score(X_test, y_test)
        # plt.scatter(y_test[:, -1], predicted_values, s=3, marker='o',label='Original')
        # plt.scatter(y_test[:, -1], predicted_values2, s=3, marker='x',label='Shuffled')
        # plt.title('R2: '+ str(R2) +'\n''R2_shuffled: '+ str(R2_shuffled) )

        R2_list.append(R2)
        R2_shuffled_list.append(R2_shuffled)
        coefficient_list.append(model.coef_)
        intercept_list.append(model.intercept_)
    return R2_list, R2_shuffled_list, coefficient_list, intercept_list