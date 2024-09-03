import pandas as pd
import numpy as np
import os

def significance(stat):
    """
    Calculate the number of stars based on the p-value.
    """
    if stat>0.05:
        stars=0
    elif 0.01<=stat<0.05:
        stars=1
    elif 0.001<=stat<0.01:
        stars=2
    elif 0.0001<=stat<0.001:
        stars=3
    elif stat<0.0001:
        stars=4
    return stars


def save_stats_data(annotations, Dir, filename):
    """
    Save the statistical annotations (testnamem p-value, n#) data to a csv file and the total data to a numpy file.
    """
    df = pd.DataFrame({'test_name': [], 'group': [], 'test_statistic': [], 'p-value': [], 'n_value': [],'stars':[]})
    for i in range(len(annotations)):
        groups = []
        numbers = []
        for j in range(len(annotations[i].structs)):
            groups.append(annotations[i].structs[j]['label'])
            numbers.append(annotations[i].structs[j]['group_data'].shape[0])
        df = pd.concat([df, pd.DataFrame(
            {'test_name': [annotations[i].data.test_description + ' ' + str(annotations[i].data.alpha)], 'group': [groups], 'test_statistic': [annotations[i].data.stat_value],
             'p-value': [annotations[i].data.pvalue], 'n_value': [numbers], 'stars':[significance(annotations[i].data.pvalue)]})], ignore_index=True)
    df.to_csv(os.path.join(Dir, filename + '_stats.csv'), index=False)
    total_data = {}
    for i in range(len(annotations)):
        for j in range(len(annotations[i].structs)):
            total_data[annotations[i].structs[j]['label']] = np.array(annotations[i].structs[j]['group_data'])
    np.save(os.path.join(Dir, filename + '_total_data.npy'), total_data, allow_pickle=True)

    return 1