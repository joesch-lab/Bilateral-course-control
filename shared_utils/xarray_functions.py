import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
import dask 
dask.config.set(**{'array.slicing.split_large_chunks': True})
from scipy.stats import spearmanr
# from scipy.signal import correlate
from statsmodels.tsa.stattools import acf, ccf
from scipy.stats import sem
from scipy.stats import entropy

from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.preprocessing import scale
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.cm as cm
from scipy.signal import savgol_filter
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def adjust_frame_rate(data, new_frame_rate=100):
    return data.resample(time='{}ms'.format(1000//new_frame_rate)).mean().interpolate_na('time')


def z_score(data):
    func = lambda x: StandardScaler().fit_transform(x.T).T
    return xr.apply_ufunc(func, data, dask='parallelized')


def stack_by_trials(data, dropna=True):
    dims = list(data.dims)
    dims_to_stack = dims.copy()
    for x in dims:
        if x in ['time']:
            dims_to_stack.remove(x)
    data=data.stack(all_trials = dims_to_stack)
    data=data.assign_coords(all_trials = np.arange(data.all_trials.size)).\
    drop_vars(dims_to_stack).\
    transpose(..., 'time')
    if dropna:
        data = data.dropna('all_trials', how='all') # drop trials with all NaNs
    return data


def sort_along_dim(data_to_sort, sort_by, dim):
    func = lambda x, y: x[np.argsort(y), :]
    for d in list(data_to_sort.dims):
        if 'time' in d:
            time_dim = d
    return xr.apply_ufunc(func, data_to_sort, sort_by,
                        input_core_dims=[[dim, time_dim], [dim]],  # list with one entry per arg
                        output_core_dims=[[dim, time_dim]],
                        output_dtypes=((float),),  
                        vectorize=True,  # loop over non-core dims
                        dask='parallelized')


def shuffle_data(a):
    rng = np.random.default_rng()
    func = lambda x: rng.permuted(x, axis=-1)
    return xr.apply_ufunc(func, a, dask='parallelized')

def shuffle_data_along_axis(a):
    def _shuffle(x):
        rng = np.random.default_rng()
        rng.shuffle(x, axis=-1)
        return x
    return xr.apply_ufunc(_shuffle, a, dask='parallelized')

def roll_data(a, shift=None):
    if shift is None:
        shift = np.random.randint(20, 40)
    func = lambda x: np.roll(x, shift, axis=-1)
    return xr.apply_ufunc(func, a, dask='parallelized')


def spearman_correlation(a, b, dim):
    func = lambda a, b: spearmanr(a, b, nan_policy='propagate').correlation
    return xr.apply_ufunc(func, a, b,
    input_core_dims=[[dim], [dim]],  # list with one entry per arg
    exclude_dims=set((dim,)),
    vectorize=True,  # loop over non-core dims
    output_dtypes = [float],
    dask='parallelized')


def autocorrelation(a, dim, output_size=120):
    func = lambda a: np.concatenate((acf(a[::-1], nlags=a.shape[-1]-1)[::-1], acf(a, nlags=a.shape[-1]-1)))
    return xr.apply_ufunc(func, a,
            input_core_dims=[[dim]],  # list with one entry per arg
            output_core_dims=[[dim]],  
            exclude_dims=set((dim,)),  # dimensions allowed to change size.  
            output_dtypes=((float),),
            vectorize=True,  # loop over non-core dims
            output_sizes={'new_time': output_size},
            dask='parallelized')


def crosscorrelation(a, b, dim, output_size=120):
    # func = lambda a: correlate(a, a, mode='same')
    func = lambda a,b: np.concatenate((ccf(a[::-1], b[::-1])[::-1], ccf(a, b)))
    return xr.apply_ufunc(func, a, b,
            input_core_dims=[[dim], [dim]],  # list with one entry per arg
            output_core_dims=[[dim]],
            exclude_dims=set((dim,)),  # dimensions allowed to change size.    
            output_dtypes=((float),),
            vectorize=True,  # loop over non-core dims
            output_sizes={'new_time': output_size},
            dask='parallelized')


def savgol_smoothening_courtship(data, window_length=31, order=2, dim='time'):
    data['angular_speed_m'] = _savgol(data.angular_speed_m, window_length=window_length, polyorder=order, dim=dim).compute()
    data['speed_m'] = _savgol(data.speed_m, window_length=window_length, polyorder=order, dim=dim).compute()
    data['angular_speed_f'] = _savgol(data.angular_speed_f, window_length=window_length, polyorder=order, dim=dim).compute()
    data['speed_f'] = _savgol(data.speed_f, window_length=window_length, polyorder=order, dim=dim).compute()
    return data


def savgol_smoothening(data, window_length=31, order=2, dim='time'):
    data['angular_speed'] = _savgol(data.angular_speed, window_length=window_length, polyorder=order, dim=dim).compute()
    data['speed'] = _savgol(data.speed, window_length=window_length, polyorder=order, dim=dim).compute()
    return data


def _savgol(data, window_length=31, polyorder=2, dim='time'):
    ### nan values are removed(for applying savgol filter) and appeneded back to the array
    func = lambda a : np.append(savgol_filter(a[~np.isnan(a)], window_length=window_length, polyorder=polyorder, axis=-1), \
                                np.empty(np.count_nonzero(np.isnan(a)))*np.nan)
    return xr.apply_ufunc(func, data,
            input_core_dims=[[dim]],  # list with one entry per arg
            output_core_dims=[[dim]],
            exclude_dims=set((dim,)),  # dimensions allowed to change size.    
            output_dtypes=((float),),
            vectorize=True,  # loop over non-core dims
            output_sizes={'time': data.time.size},
            dask='parallelized')

def xr_sem(a, dim):
    dims = list(a.dims)
    for d in dims:
        if 'time' in d:
            time_dim = d
    a = a.transpose(..., time_dim)
    func = lambda a: sem(a, nan_policy='omit', axis=0)
    return xr.apply_ufunc(func, a,
                exclude_dims=set((dim,)),
                input_core_dims=[[dim, time_dim]],  # list with one entry per arg
                output_core_dims=[[time_dim]], 
                output_dtypes=((float),),
                vectorize=True,  # loop over non-core dims
                dask='parallelized',
                dask_gufunc_kwargs={'allow_rechunk': True})


def univeral_log(data):
    func = lambda x: np.log2(x) if x>1 else (-np.log2(-x) if x<-1 else 0)
    return xr.apply_ufunc(func, data,
            vectorize=True,  # loop over non-core dims
            dask='parallelized',
            dask_gufunc_kwargs={'allow_rechunk': True})


def create_new_direction_dim(data, stim_start=720):
    direction_data = data.\
                    isel(time=slice(stim_start, None)).\
                    pinwheel_angle.diff('time').isel(time=0).compute() # gives the direction of the pinwheel motion
    print(np.unique(direction_data.values))
    cw_data = data.where(direction_data>0)
    ccw_data = data.where(direction_data<0)
    newdata = xr.concat([cw_data, ccw_data], dim='direction')
    newdata = newdata.assign_coords(direction=[-1, 1])
    return newdata


def change_for_expected_direction(data):
    """
    This function changes the direction of the angular speed data to match the expected direction of the angular speed.
    data: xarray dataset
    variable: variable name to change the direction of
    """
    midpoint = data.time.shape[0]//2
    mid = data.time.values[midpoint]
    data['orientation'] = data['orientation'] - data['orientation'].loc[dict(time=mid)]
    data['angular_speed'].loc[dict(direction=1)] = data['angular_speed'].loc[dict(direction=1)] * -1
    data['angular_speed'].loc[dict(direction=-1)] = data['angular_speed'].loc[dict(direction=-1)]
    data['orientation'].loc[dict(direction=1)] = data['orientation'].loc[dict(direction=1)] * -1
    data['orientation'].loc[dict(direction=-1)] = data['orientation'].loc[dict(direction=-1)]
    try:
        data['saccades'].loc[dict(direction=1)] = data['saccades'].loc[dict(direction=1)] * -1
        data['saccades'].loc[dict(direction=-1)] = data['saccades'].loc[dict(direction=-1)]
    except:
        print('No saccades')
    return data


def change_for_expected_direction_multiindex(data):
    """
    This function changes the direction of the angular speed data to match the expected direction of the angular speed.
    data: xarray dataset
    variable: variable name to change the direction of
    """
    midpoint = data.time.shape[0]//2
    mid = data.time.values[midpoint]

    for i in range(data.motion.direction.shape[0]):
        data['orientation'].loc[dict(motion=data.motion[i])] = data['orientation'].loc[dict(motion=data.motion[i])] - data['orientation'].loc[dict(motion=data.motion[i], time=mid)]
        if data.motion.direction[i] == 1:
            data['angular_speed'].loc[dict(motion=data.motion[i])] = data['angular_speed'].loc[dict(motion=data.motion[i])] * -1
            data['orientation'].loc[dict(motion=data.motion[i])] = data['orientation'].loc[dict(motion=data.motion[i])] * -1
        elif data.motion.direction[i] == -1:
            data['angular_speed'].loc[dict(motion=data.motion[i])] = data['angular_speed'].loc[dict(motion=data.motion[i])]
            data['orientation'].loc[dict(motion=data.motion[i])] = data['orientation'].loc[dict(motion=data.motion[i])]

    return data


def change_for_expected_direction_courtship(data):
    """
    This function changes the direction of the angular speed data to match the expected direction of the angular speed.
    data: xarray dataset
    variable: variable name to change the direction of
    """
    midpoint = data.time.shape[0]//2
    mid = data.time.values[midpoint]
    print(mid)
    data['ori_smooth_m'] = data['ori_smooth_m'] - data['ori_smooth_m'].loc[dict(time=mid)]
    data['angular_speed_m'].loc[dict(direction=1)] = data['angular_speed_m'].loc[dict(direction=1)]
    data['angular_speed_m'].loc[dict(direction=-1)] = data['angular_speed_m'].loc[dict(direction=-1)] * -1
    data['ori_smooth_m'].loc[dict(direction=1)] = data['ori_smooth_m'].loc[dict(direction=1)]
    data['ori_smooth_m'].loc[dict(direction=-1)] = data['ori_smooth_m'].loc[dict(direction=-1)] * -1

    data['ori_smooth_f'] = data['ori_smooth_f'] - data['ori_smooth_f'].loc[dict(time=mid)]
    data['angular_speed_f'].loc[dict(direction=1)] = data['angular_speed_f'].loc[dict(direction=1)]
    data['angular_speed_f'].loc[dict(direction=-1)] = data['angular_speed_f'].loc[dict(direction=-1)] * -1
    data['ori_smooth_f'].loc[dict(direction=1)] = data['ori_smooth_f'].loc[dict(direction=1)]
    data['ori_smooth_f'].loc[dict(direction=-1)] = data['ori_smooth_f'].loc[dict(direction=-1)] * -1
    return data


def get_position_histogram(data_x, data_y,dim, bins=100, range=[0,300]):
    func = lambda x, y: np.histogram2d(x, y, bins=bins, range=[range, range], density=True)[0]
    return xr.apply_ufunc(func, data_x, data_y,
                            input_core_dims=[[dim], [dim]],  # list with one entry per arg
                            exclude_dims=set((dim,)),
                            output_core_dims=[['x', 'y']],
                            vectorize=True,  # loop over non-core dims
                            output_dtypes = ((float), (float)),
                            output_sizes={'x': bins, 'y': bins},
                            dask='parallelized', 
                            dask_gufunc_kwargs={'allow_rechunk':True})


def get_histogram(data_ori, dim, bins=180, range=[0,180]):
    func = lambda x: np.histogram(x, bins=bins, range=range, density=True)[0]
    return xr.apply_ufunc(func, data_ori,
            input_core_dims=[[dim]],  # list with one entry per arg
            exclude_dims=set((dim,)),
            output_core_dims=[['angles']],
            vectorize=True,  # loop over non-core dims
            output_dtypes = ((float)),
            output_sizes={'angles': bins},
            dask='parallelized', 
            dask_gufunc_kwargs={'allow_rechunk':True})


def pseudo_MI_2d(a, b, c, d, dim, shift=5):
    func = lambda a,b,c,d: entropy(get_position_histogram(a, b, dim, 200, [0,1000]).flatten(), base=2) - entropy(get_position_histogram(c, d, dim, 200, [0,1000]).flatten(), base=2)
    
    return xr.apply_ufunc(func, a, b, c.shift(time=shift), d.shift(time=shift),
        input_core_dims=[[dim],[dim],[dim],[dim]],  # list with one entry per arg
        exclude_dims=set((dim,)),
        vectorize=True,  # loop over non-core dims
        output_dtypes = [float],
        dask='parallelized', 
        dask_gufunc_kwargs={'allow_rechunk':True})


def male_to_female_analysis(pos_x_m, pos_y_m, pos_x_f, pos_y_f, ori_360_m, ori_360_f):
    """compute relative distance and angles, to plot the distributions as evidence of tracking"""
    male_to_female_x = np.subtract(pos_x_f, pos_x_m)
    male_to_female_y = np.subtract(pos_y_f, pos_y_m)
    male_to_female_dist = np.sqrt(np.add(np.square(male_to_female_x), np.square(male_to_female_y)))
    # the angle of the line joining the male center to the female center
    male_to_female_dir = (np.multiply(np.arctan2(male_to_female_y*-1, male_to_female_x) / np.pi , 180) + 90) % 360
    # the orientation of the male wrt the line joining the centers of male and female
    male_ori_wrt_female = np.subtract(male_to_female_dir, ori_360_m) % 360
    male_ori_wrt_female = 180 + np.abs(np.subtract(male_ori_wrt_female, 180))
    # the orientation of the female wrt the line joining the centers of male and female
    female_ori_wrt_male = np.add(np.subtract(male_to_female_dir, ori_360_f), 180) % 360
    female_ori_wrt_male = 180 + np.abs(np.subtract(female_ori_wrt_male, 180))
    return male_to_female_dist, male_ori_wrt_female, female_ori_wrt_male


def male_to_female_analysis_xr(data):
    func = male_to_female_analysis(data.pos_x_m, data.pos_y_m, data.pos_x_f, data.pos_y_f, data.ori_360_m, data.ori_360_f)
    return xr.apply_ufunc(func, data,
        input_core_dims=[['time']],  # list with one entry per arg
        output_core_dims=[['time']],
        output_dtypes=((float),),
        vectorize=True,  # loop over non-core dims
        dask='parallelized')


# sort a list of strings based on the letters in the string
def sort_list_of_strings(list_of_strings):
    return sorted(list_of_strings, key=lambda v: v.upper())


def _get_mean_saccade(angular_speed, saccade_data):
    syn_saccades = []
    anti_saccades = []
    ind = np.argwhere(saccade_data>0)
    l = angular_speed.shape[1]
    for i, j in ind:
        if 6 < j < l-6:
            syn_saccades.append(angular_speed[i, j-5:j+5])
    ind = np.argwhere(saccade_data<0)
    for i, j in ind:
        if 6 < j < l-6:
            anti_saccades.append(angular_speed[i,j-5:j+5])
    output = np.empty((2, 10))*np.nan
    if len(syn_saccades)>0:
        output[0,:] = np.mean(syn_saccades, axis=0)
    if len(anti_saccades)>0:
        output[1,:] = np.mean(anti_saccades, axis=0)
    return output


def get_mean_saccade(data1, data2):
    return xr.apply_ufunc(_get_mean_saccade, data1, data2,
                        input_core_dims=[['trials', 'time'], ['trials', 'time']],  # list with one entry per arg
                        output_core_dims=[['x', 'y']],
                        output_dtypes=((float), (float)),
                        output_sizes={'x': 2, 'y':10},  
                        vectorize=True,  # loop over non-core dims
                        dask='parallelized')


# return with cluster labels
def _cluster(data, gmm, pca, n_samples=100):
    good_indices = np.where(~np.isnan(data).all(axis=1))[0]
    final_data = data[~np.any(np.isnan(data), axis=1), :]
    scaled_data = scale(final_data.T).T

    # scaled_data = pca.transform(scaled_data)
    # dims_to_keep = 10
    # scaled_data = scaled_data[:, :dims_to_keep]

    result = gmm.predict(scaled_data)
    cluster_ids = np.empty(data.shape[0]) * np.nan
    cluster_ids[good_indices] = result
    output = np.empty((gmm.means_.shape[0], n_samples, scaled_data.shape[1])) * np.nan
    for i in range(gmm.means_.shape[0]):
        x = final_data[np.argwhere(result==i).flatten(), :]
        output[i, :x.shape[0], :] = x
    return cluster_ids, output


def plot_silhouette_values(X, cluster_labels, n_clusters=3, name=''):
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =",n_clusters,"The average silhouette_score is :",silhouette_avg,)
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    y_lower = 10
    fig, ax = plt.subplots()
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper), 0, ith_cluster_silhouette_values, facecolor=color, edgecolor=color, alpha=0.7,)
        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")
    
    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xlim([-0.5, 1])
    name = name + str(n_clusters)
    ax.set_title("Avg silhouette score for {}: {}".format(name, str(silhouette_avg)))
    fig.savefig(name+'_silhouette.png')

    return silhouette_avg


def cluster_tracking_data(data, n_clusters=4, feature='pinwheel_contrast', trial_dim='blocked_trials', gmm=None, pca=None):
    final_data = data.stack(x=(feature, trial_dim)).transpose(..., 'new_time').compute()
    final_data = scale(final_data.T).T
    if not gmm:
        # pca = PCA(n_components=10)
        # final_data = pca.fit_transform(final_data)
        # dims_to_keep = 5
        # final_data = final_data[:, :dims_to_keep]
        gmm = GaussianMixture(n_components=n_clusters, random_state=2)
        # gmm = KMeans(n_clusters=n_clusters)
        final_data = final_data[~np.any(np.isnan(final_data), axis=1), :]
        cluster_labels = gmm.fit_predict(final_data)
        plot_silhouette_values(final_data, cluster_labels, n_clusters=n_clusters)
    return xr.apply_ufunc(_cluster, data, gmm, pca, data[trial_dim].shape[0],
                        input_core_dims=[[trial_dim, 'new_time'], [], [], []],  # list with one entry per arg
                        exclude_dims=set((trial_dim,)),
                        output_core_dims=[['x'],['clusters', 'x', 'new_time']],
                        output_dtypes=[np.float64, np.float64],
                        output_sizes={'clusters': n_clusters, 'x':data[trial_dim].shape[0]},  
                        vectorize=True,  # loop over non-core dims
                        dask='parallelized',
                        dask_gufunc_kwargs={'allow_rechunk':True})


def _remove_saccades(angular_speed, orientation, saccades, remove='all', stim_period=300, window=[5,5], filler_type='zero'):
    new_angular_speed = angular_speed.copy()
    new_orientation = orientation.copy()
    for i in range(orientation.shape[0]):
        if remove == 'anti' or remove == 'all':
            anti_indices = np.argwhere(saccades[i, stim_period:] < 0)
            for index in anti_indices:
                if index[0] > window[0] and index[0] < (orientation.shape[1]-stim_period) - window[1]:
                    new_orientation[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = np.zeros(window[0] + window[1])
                    if filler_type == 'zero':
                        filler = np.zeros(window[0] + window[1])
                    elif filler_type == 'linear':
                        filler = np.interp(np.arange(0, window[0] + window[1], 1), [-1, window[0] + window[1] + 1],
                                            [new_angular_speed[i, index[0] + stim_period - window[0] - 1],
                                            new_angular_speed[i, index[0] + stim_period - window[0] + 1]])
                        filler = np.random.rand(window[0] + window[1]) * 25 + filler
                    new_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = filler
        if remove == 'syn' or remove == 'all':
            anti_indices = np.argwhere(saccades[i, stim_period:] > 0)
            for index in anti_indices:
                if index[0] > window[0] and index[0] < (orientation.shape[1]-stim_period) - window[1]:
                    new_orientation[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = np.zeros(window[0] + window[1])
                    if filler_type == 'zero':
                        filler = np.zeros(window[0] + window[1])
                    elif filler_type == 'linear':
                        filler = np.interp(np.arange(0, window[0] + window[1], 1), [-1, window[0] + window[1] + 1],
                                            [new_angular_speed[i, index[0] + stim_period - window[0] - 1],
                                            new_angular_speed[i, index[0] + stim_period - window[0] + 1]])
                        filler = np.random.rand(window[0] + window[1]) * 25 + filler
                    new_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = filler
    return new_orientation, new_angular_speed


def _keep_saccades(angular_speed, orientation, saccades, keep='all', stim_period=300, window=[5,5]):
    new_angular_speed = np.zeros(angular_speed.shape)
    new_orientation = np.zeros(orientation.shape)
    for i in range(orientation.shape[0]):
        if keep == 'anti' or keep == 'all':
            anti_indices = np.argwhere(saccades[i, stim_period:] < 0)
            for index in anti_indices:
                if index[0] > window[0] and index[0] < (orientation.shape[1]-stim_period) - window[1]:
                    new_orientation[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = \
                    orientation[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]]
                    new_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = \
                    angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]]
        if keep == 'syn' or keep == 'all':
            anti_indices = np.argwhere(saccades[i, stim_period:] > 0)
            for index in anti_indices:
                if index[0] > window[0] and index[0] < (orientation.shape[1]-stim_period) - window[1]:
                    new_orientation[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]] = \
                    orientation[i, index[0]+stim_period-window[0]:index[0]+stim_period+window[1]]
                    new_angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]] = \
                    angular_speed[i, index[0] + stim_period - window[0]:index[0] + stim_period + window[1]]
    return new_orientation, new_angular_speed


def remove_saccades_xarray(angular_speed, orientation, saccades, remove='all'):
    return xr.apply_ufunc(_remove_saccades, angular_speed, orientation, saccades,
                    input_core_dims=[['trials', 'time'], ['trials', 'time'], ['trials', 'time']],  # list with one entry per arg
                    output_core_dims=[['trials', 'time'], ['trials', 'time']],
                    kwargs={'remove': remove},
                    output_dtypes=((float), (float)), 
                    vectorize=True,  # loop over non-core dims
                    dask='parallelize')


def keep_saccades_xarray(angular_speed, orientation, saccades, keep='all'):
    return xr.apply_ufunc(_remove_saccades, angular_speed, orientation, saccades,
                    input_core_dims=[['trials', 'time'], ['trials', 'time'], ['trials', 'time']],  # list with one entry per arg
                    output_core_dims=[['trials', 'time'], ['trials', 'time']],
                    output_dtypes=((float), (float)), 
                    kwargs={'keep': keep},
                    vectorize=True,  # loop over non-core dims
                    dask='parallelize')
