import copy
import os
import json
import pandas as pd
import numpy as np
import pandas as pd
import pickle
from datetime import datetime

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class fly_exp_data:
    def __init__(self,json_name,json_data,csv_name,csv_data,csv_stimulus_name,csv_stimulus_data):
        self.json_name = json_name
        self.csv_name = csv_name
        self.track = csv_data
        self.fly = json_data
        self.stimulus = csv_stimulus_data
        self.stimulus_name = csv_stimulus_name

def load_json(Dir,find='',ignore='zzzz'):
    for files in os.listdir(Dir):
        if files.endswith('json') and find in files and ignore not in files:
            path = os.path.join(Dir, files)
            with open(path, 'r') as f:
                try:
                    data = json.load(f)
                except:
                    p2jfile_hemi(Dir)
                    f = open(path, 'r')
                    data = json.load(f)
    return files, data

def load_csv(Dir,find='',ignore='', fix = 'N'):
    file=[]
    data_comb=[]
    for files in os.listdir(Dir):
        if files.endswith('csv') and find in files and ignore not in files and 'stimulus' not in files and 'old' not in files and 'notfixed' not in files:
            print(files)
            path = os.path.join(Dir, files)
            data = pd.read_csv(path, on_bad_lines='skip')
            data = data.dropna(how='all')
            if fix =='Y':
                print('fixing')
                data = fix_csv(data)
            data_comb.append(data)
            file.append(files)
    return file, data_comb

def load_stimulus_csv(Dir,find='',ignore='', fix = 'N'):
    file=[]
    data_comb=[]
    for files in os.listdir(Dir):
        if files.endswith('csv') and find in files and ignore not in files and 'stimulus_conditions' in files and 'old' not in files and 'notfixed' not in files:
            path = os.path.join(Dir, files)
            data = pd.read_csv(path, on_bad_lines='skip', skip_blank_lines=True, header=0)
            data_comb.append(data)
            file.append(files)
    return file, data_comb

def search_data(Dir,data_loaded,find='',ignore='zzzz', fix_csv='N'):
    print(Dir, 'here')
    """
    str(name of directory) -> list(fly_exp_data object)
    assumption : each folder has data for one fly and hence one .json file and one .p file
                there might be multiple .csv files,
                each holding data for one experiment instance
    output : each element of a list is one such folder data - one .json file (organised as a dict object)
            and one(or more) .csv file (as a 2D pandas dataframe)
    """
    for files in os.listdir(Dir):
        path = os.path.join(Dir,files)
        if os.path.isdir(path):
            search_data(path,data_loaded)
        else:
            if files.endswith('csv') and find in files and ignore not in files and 'stimulus' not in files and 'old' not in files and 'notfixed' not in files:
                print('data file: ', files)
                json_name, json_data = load_json(Dir,ignore=ignore)
                # if json_data == {}:
                #     print('Error')
                #     convert_pickle_to_json(files.strip('.json')+'.p')
                #     json_name, json_data = load_json(Dir, find=find, ignore=ignore)
                csv_name, csv_data = load_csv(Dir,find=find,ignore=ignore, fix=fix_csv)
                csv_stimulus_name, csv_stimulus_data = load_stimulus_csv(Dir,find=find,ignore=ignore, fix=fix_csv)
                fly = fly_exp_data(json_name,json_data,csv_name,csv_data, csv_stimulus_name, csv_stimulus_data)
                print(csv_name)
                data_loaded.append(fly)
    return data_loaded

def search_old_data(Dir,data_loaded,find='',ignore='zzzz'):
    #print(Dir)
    """
    str(name of directory) -> list(fly_exp_data object)
    assumption : each folder has data for one fly and hence one .json file and one .p file
                there might be multiple .csv files,
                each holding data for one experiment instance
    output : each element of a list is one such folder data - one .json file (organised as a dict object)
            and one(or more) .csv file (as a 2D pandas dataframe)
    """
    for files in os.listdir(Dir):
        path = os.path.join(Dir, files)
        if os.path.isdir(path):
            search_old_data(path, data_loaded)
        else:
            if files.endswith('.csv') and find in files and ignore not in files:
                json_name = 'blah'
                json_data = {'name':'blah'}
                csv_name, csv_data = load_csv(Dir, find=find, ignore=ignore)
                fly = fly_exp_data(json_name, json_data, csv_name, csv_data)
                data_loaded.append(fly)
    return data_loaded

def p2jfile(Dir):
    # Dir = r'D:\Roshan\Experiments\Optogenetics\P1\P1_Courtship\160920\P1_csChrimson_size_8'
    for files in os.listdir(Dir):
        if files.endswith('.p'):
            path = os.path.join(Dir, files)
            data = convert_pickle_to_json(path)
            # for stuff in data['exp_params'][0]['stim_params']:
            #     stuff['stim'][2] = list(stuff['stim'][2])
            filename_json = path.replace('.p', '.json')
            json_file = open(filename_json, 'w')
            json.dump(data, json_file, cls=NumpyEncoder)
            json_file.close()
    return 1

def p2jfile_new(Dir):
    # Dir = r'D:\Roshan\Experiments\Optogenetics\P1\P1_Courtship\160920\P1_csChrimson_size_8'
    for files in os.listdir(Dir):
        if files.endswith('.p'):
            path = os.path.join(Dir, files)
            data = convert_pickle_to_json(path)
            # for stuff in data['exp_params'][0]['stim_params']:
            #     stuff['stim'][2] = list(stuff['stim'][2])
            filename_json = path.replace('.p', '.json')
            json_file = open(filename_json, 'w')
            json.dump(data, json_file, cls=NumpyEncoder)
            json_file.close()
    return 1

def convert_pickle_to_json(filename):
    '''
    str(filename(p-file)) -> json dict
    the ourput can be directly stored as a json file by using json.dump(read_pickle(filename))
    :param filename: str
    :return: list off dicts(json formate data)
    '''
    datastore = []
    f = open(filename,'rb')
    #pickle.seek(0) ## sets pointer at the beginning of the file
    while(True):
        try:
            data = pickle.load(f)
            datastore.append(data)
        except Exception as error:
            print('Read the file')
            break
    
    if 'id' in datastore[0]:
        datadict = {**datastore[0], **{'exp_params':[]}}
    else:
        id = os.path.basename(filename).split('_')[-1].split('.')[0]
        sex='M'
        age=2
        strain = os.path.basename(filename).split('_')[0]
        datadict = {**{'id': id, 'sex': sex, 'age': age, 'strain': strain}, **{'exp_params': []}}
        t = datetime.fromtimestamp(os.path.getctime(filename))
        datadict['exp_params'].append({**{'res': (720, 720), 'fps': 60, 'time': [t.time().hour, t.time().minute, t.time().second, t.time().microsecond], 
                                          'date': [t.year, t.month, t.day], 
                                          'bg': 1, 'fg': -1, 'src_file': filename},**{'stim_params':[]}})
    exp_number = 0
    for dicts in datastore:
        if 'src_file' in dicts:
            datadict['exp_params'].append({**dicts, **{'stim_params':[]}})
            exp_number += 1
            continue
        elif 'on_frames' in dicts:
            print(datadict['exp_params'])
            datadict['exp_params'][exp_number - 1]['stim_params'].append(dicts)
        else:
            continue
    return datadict

def fix_csv(a):
    pos_x = copy.deepcopy(a.index)
    pos_y = copy.deepcopy(a['pos_x'])
    ori = copy.deepcopy(a['pos_y'])
    timestamp = copy.deepcopy(a['ori'])
    direction = copy.deepcopy(a['timestamp'])
    a = a.drop(columns=['direction'])
    a['pos_x'] = pos_x
    a['pos_y'] = pos_y
    a['ori'] = ori
    a['timestamp'] = timestamp
    a['direction'] = direction
    a.reset_index(drop=True, inplace=True)
    a.to_csv()
    return a

def p2jfile_hemi(Dir):
    # Dir = r'D:\Roshan\Experiments\Optogenetics\P1\P1_Courtship\160920\P1_csChrimson_size_8'
    for files in os.listdir(Dir):
        if files.endswith('.p'):
            path = os.path.join(Dir, files)
            data = convert_pickle_to_json_hemi(path)
            # for stuff in data['exp_params'][0]['stim_params']:
            #     stuff['stim'][2] = list(stuff['stim'][2])
            filename_json = path.replace('.p', '.json')
            json_file = open(filename_json, 'w')
            json.dump(data, json_file)
            json_file.close()
    return 1

def convert_pickle_to_json_hemi(filename):
    '''
    str(filename(p-file)) -> json dict
    the ourput can be directly stored as a json file by using json.dump(read_pickle(filename))
    :param filename: str
    :return: list off dicts(json formate data)
    '''
    datastore = []
    f = open(filename, 'rb')
    #pickle.seek(0) ## sets pointer at the beginning of the file
    while(True):
        try:
            data = pickle.load(f)
            datastore.append(data)
        except:
            print('Read the file')
            break

    datadict = {**datastore[0], **{'exp_params': []}}
    exp_number = 0
    i = 0

    print(len(datastore))
    while i < (len(datastore)//2)*2:
        if 'src_file' in datastore[i]:
            datadict['exp_params'].append({**datastore[i], **{'stim_params': []}})
            exp_number += 1
        elif 'stim_size' in datastore[i]:
            datadict['exp_params'][exp_number - 1]['stim_params'].append([datastore[i], datastore[i + 1]])
            i += 1
        else:
            pass
        i += 1

    return datadict

def p2jfile_latest(Dir):
    for file in os.listdir(Dir):
        if file.endswith('.p'):
            break
    filename = os.path.join(Dir, file)
    f = open(filename, 'rb')
    newdata = pickle.load(f)
    while (True):
        try:
            newdata['exp_params'][0]['stim_params'].append(pickle.load(f))
        except:
            print('Read the file')
            break

    newfilename = os.path.join(Dir, file[:-2] + '.json')
    f = open(newfilename, 'w')
    json.dump(newdata, f, cls=NumpyEncoder)
    f.close()
    return 1

def stim_params_list_json(flydata, param, index):
    """
    Get the list of values for a given parameter from flydata
    :param flydata: object of class FlyData
    :type flydata: FlyData
    :param param: parameter name
    :type param: str
    :param index: index of the experiment, default is 0
    :type index: int
    :return: list of values for the given parameter
    :rtype: list
    """
    if type(flydata.fly['exp_params'][index]['stim_params'][0]) == dict:
        param_list = []
        for i in range(len(flydata.fly['exp_params'][index]['stim_params'])):
            param_list.append(flydata.fly['exp_params'][index]['stim_params'][i][param])
    elif type(flydata.fly['exp_params'][index]['stim_params'][0]) == list:
        param_list = []
        for j in range(len(flydata.fly['exp_params'][index]['stim_params'][0])):
            params = []
            for i in range(len(flydata.fly['exp_params'][index]['stim_params'])):
                params.append(flydata.fly['exp_params'][index]['stim_params'][i][j][param])
            param_list.append(params)
    return param_list

if __name__ == "__main__":
    try:
        search_data(r"D:/Roshan/Zzzz")
    except:
        print("No such directory")