import numpy as np

def subplot_arrangement(n):
    """
    Calculate the number of rows and columns for a subplot arrangement.
    :param n: The number of subplots.
    :type n: int
    :returns: The number of rows and columns.
    :rtype: (int, int)
    """
    if n<5:
        return 1, n
    else:
        m = np.sqrt(n)
        i = int(m)
        j = int(round(m))
        if i*j == n:
            return i,j
        elif i*j > n:
            return i, j-1
        elif i*j < n:
            return i+1, j
    return 1


def make_labels(conditions):
    """
    create labels to be used in plotting, from the stimulus conditions.
    :param conditions: The stimulus conditions
    :type conditions: list of str or list of list of str
    :returns: The labels
    :rtype: list of str
    """
    names = []
    if len(conditions[0])==1:
        names = [str(x) for x in list(conditions[0][0].astype(int))]
    else:
        for i in range(len(conditions)):
            new_name = str(conditions[i][0])
            for j in range(1, len(conditions[i])):
                new_name = new_name + '_' + str(conditions[i][j])
            names.append(new_name)
    return names
