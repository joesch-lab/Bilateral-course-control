import numpy as np
from functools import reduce


def indices_by_condition(conditions, ignore_values=None):
    """
    Take an array of conditions and return the indices of the trials that belong to each condition.
    :param conditions: The stimulus conditions or experiment conditions.
    :type conditions: list of numpy.ndarray
    :returns: The indices of the trials that belong to each condition and the sorted unique conditions.
    :rtype: (list of numpy.ndarray, list)
    """
    if not ignore_values:
        ignore_values = [[None]] * len(conditions)
    index = []
    sorted_conditions = []
    for i in range(len(conditions)):
        unique_values = list(np.unique(conditions[i]))
        unique_values = [x for x in unique_values if x not in ignore_values[i]]
        c1 = np.sort(unique_values)
        index2 = []
        for val in c1:
            index3 = np.argwhere(conditions[i]==val).flatten()
            index2.append(index3)
        index.append(index2)
        sorted_conditions.append(c1)

    if len(index) == 1:
        return index[0], [sorted_conditions]
    elif len(index) == 2:
        indices = []
        condns = []
        for i in range(len(index[0])):
            for j in range(len(index[1])):
                common_indices = np.intersect1d(index[0][i], index[1][j])
                if common_indices.size != 0:
                    indices.append(common_indices)
                    condns.append([sorted_conditions[0][i], sorted_conditions[1][j]])
        return indices, condns
    elif len(index) == 3:
        indices = []
        condns = []
        for i in range(len(index[0])):
            for j in range(len(index[1])):
                for k in range(len(index[2])):
                    common_indices = reduce(np.intersect1d, (index[0][i], index[1][j], index[2][k]))
                    if common_indices.size != 0:
                        indices.append(common_indices)
                        condns.append([sorted_conditions[0][i], sorted_conditions[1][j], sorted_conditions[2][k]])
        return indices, condns
    elif len(index) == 4:
        indices = []
        condns = []
        for i in range(len(index[0])):
            for j in range(len(index[1])):
                for k in range(len(index[2])):
                    for l in range(len(index[3])):
                        common_indices = reduce(np.intersect1d, (index[0][i], index[1][j], index[2][k], index[3][l]))
                        if common_indices.size != 0:
                            indices.append(common_indices)
                            condns.append([sorted_conditions[0][i], sorted_conditions[1][j], sorted_conditions[2][k], sorted_conditions[3][l]])
        return indices, condns


def split_data_by_conditions(data, conditions, ignore_values=None):
    """
    Split the provided data into a list of multiple arrays based on the conditions
    :param data: List of data to be split
    :type data: list of length k (k=#data types i.e. speed, angular speed etc.) of numpy.ndarray of shape (n, m) where n is the number of trials and m is the number of frames.
    :param conditions: The stimulus conditions or experiment conditions.
    :type conditions: list of l (l=#stimulus conditions) of 1-d numpy.ndarray of shape (n,) where n is the number of trials
    :returns: The data split by conditions
    :rtype: list (length k) of list (lenght l) of numpy.ndarray"""
    indices, conditions = indices_by_condition(conditions, ignore_values=ignore_values)
    all_data = []
    for i in range(len(data)):
        new_data = []
        for j in range(len(indices)):
            new_data.append(data[i][indices[j]])
        all_data.append(new_data)
    return all_data, conditions, indices
