import numpy as np
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
from sklearn.metrics import silhouette_score


def find_k_cluster(data, num=10, k_iter = 10, iter=10, normalize="N"):
    ## num = number of clusters to aim for
    ## k_iter = number of iterations of trying to find the errors
    ## iter = maximum iterations for the algorithm for a single run
    final_data = data

    if normalize=='Y':
        final_data = scale(final_data.T).T

    k = 0
    all = []
    all_s = []

    while k < k_iter:
        j = 2
        total_distortion = []
        train_data, test_data = train_test_split(final_data, test_size = 0.25)
        total_s = []
        while j <= num:
            kmeans = KMeans(n_clusters=j, random_state=0).fit(train_data)
            result = kmeans.predict(test_data)
            total_distortion.append(kmeans.inertia_)
            #####
            total_s.append(silhouette_score(test_data, result))
            j = j + 1
        k = k + 1
        all.append(total_distortion)
        all_s.append(total_s)
    return all, all_s


def find_k_gmm_cluster(data, num=10, k_iter = 10, iter=10, normalize='N'):
    ## num = number of clusters to aim for
    ## k_iter = number of iterations of trying to find the errors
    ## iter = maximum iterations for the algorithm for a single run
    final_data = data
    if normalize=='Y':
        final_data = scale(final_data.T).T
    k = 0
    all_s1 = []
    all_s2 = []

    while k < k_iter:
        j = 2
        train_data, test_data = train_test_split(final_data, test_size = 0.25)
        total_s1 = []
        total_s2 = []
        while j <= num:
            gmm = GaussianMixture(n_components=j, random_state=0).fit(train_data)
            result = gmm.predict(test_data)
            total_s1.append(gmm.bic(test_data))
            total_s2.append(gmm.score(test_data))
            j = j + 1
        k = k + 1
        all_s1.append(total_s1)
        all_s2.append(total_s2)
    return all_s1, all_s2


def cluster_and_model_gmm(opto_response, n_clusters, normalize = 'Y'):
    final_data = np.array(opto_response)
    if normalize=='Y':
        final_data = scale(final_data.T).T
    n_clusters = n_clusters
    gmm = GaussianMixture(n_components=n_clusters, random_state=0)
    result = gmm.fit_predict(final_data)
    return gmm, result


def cluster_and_sort(data, k, minit=None, iter=10, normalize='N', kmeans=None):
    # final_data = data[:, 3+121-30+4:3+121+60+4]
    final_data = data
    if minit is None:
        minit = 'random'
    ## here I am only using the stmulus time to do the clustering
    if normalize=='Y':
        final_data = scale(final_data.T).T
    if kmeans == None:
        kmeans = KMeans(n_clusters=k)
        result = kmeans.fit_predict(final_data).reshape(-1, 1)
    else:
        result = kmeans.predict(final_data).reshape(-1, 1)
    result_new = np.vstack((result.T, data.T)).T
    order = np.argsort(result.flatten())
    result_new = result_new[order]

    values = []
    for i in range(len(result_new[:, 0])-1):
        if result_new[i, 0] != result_new[i+1, 0]:
            values.append(i+1)
    values.insert(0, 0)
    values.append(result_new.shape[0])
    mean_clusters = []
    for i in range(len(values) - 1):
        mean_clusters.append(np.array(result_new[values[i]:values[i + 1], 5:]))
    # result : the list of index of the cluster each sample belongs to
    # result_new : array of samples arranged by their respective clusters
    # mean_clusters : mean of the samples in each cluster, similar to cluster centers but not exactly
    # kmeans.cluster_centers_ : centroid of each cluster
    return result.flatten(), result_new, mean_clusters, kmeans.cluster_centers_


def cluster_and_sort_gmm(data, k, minit=None, iter=10, normalize='N', gmm = None):
    # final_data = data[:, 3+121-30+4:3+121+60+4]
    final_data = data
    if minit is None:
        minit = 'random'
    ## here I am only using the stimulus time to do the clustering
    if normalize=='Y':
        final_data = scale(final_data.T).T
    if gmm == None:
        gmm = GaussianMixture(n_components=k, random_state=0)
        result = gmm.fit_predict(final_data).reshape(-1, 1)
    else:
        result = gmm.predict(final_data).reshape(-1, 1)
    result_new = np.vstack((result.T, final_data.T)).T
    order = np.argsort(result.flatten())
    result_new = result_new[order]

    values = []
    for i in range(len(result_new[:, 0])-1):
        if result_new[i, 0] != result_new[i+1, 0]:
            values.append(i+1)
    values.insert(0, 0)
    values.append(result_new.shape[0])
    mean_clusters = []
    for i in range(len(values) - 1):
        # mean_clusters.append(np.array(result_new[values[i]:values[i + 1], 5:]))
        mean_clusters.append(np.array(result_new[values[i]:values[i + 1], 1:]))
    return result.flatten(), result_new, mean_clusters, gmm.means_
