import os
import numpy as np
import pandas as pd


def analyse_folders(Dir, filename='example', PARAMETERS=['frequency'], p=3,
                    condn_indices=[2, 1, 0], ignore_values=[], variables=[], _combine=True, _analysis=True, _saccade=True):
    result = 1
    if not os.path.exists(os.path.dirname(Dir)):
        Dir = genotype_folders_basic[Dir]
    if _combine:
        combine_data(Dir, filename=filename, sex='', PARAMETERS=PARAMETERS, saccade='Y', normal='Y')
    if _analysis:
        analysis(Dir, filename=filename, sex='', p=p, condn_indices=condn_indices, ignore_values=ignore_values, variables=variables)
    if _saccade:
        analyse_saccade(Dir, filename=filename, sex='', p=p, condn_indices=condn_indices, ignore_values=ignore_values, variables=variables)
    return result


def combine_data(Dir, filename=None, sex='', PARAMETERS = ['frequency'], saccade='Y', normal='Y'):
    if filename==None:
        filename = Dir.split(r'\'')[0].split('\\')[-1]
    total_response_data = []
    total_orientation_data = []
    total_response_speed_data = []
    total_response_angular_speed_data = []
    total_response_fwd_vel = []
    total_response_side_vel = []
    saccade_peak_list = []
    saccade_turn_list = []
    pos_list = []
    fly_id = 0
    fly_name_dict = {'filename':[], 'flyid':[]}
    flyname_file = os.path.join(Dir, 'filename.csv')

    Dir_list = glob.glob(os.path.join(Dir, '**', '*.avi'), recursive=True)
    Dir_list = [os.path.dirname(x) for x in Dir_list]
    for path in Dir_list:
        if os.path.isdir(path):
            if sex == '':
                print('here')
                data = main(path, PARAMETERS = PARAMETERS, fly_id = fly_id, saccade=saccade, normal=normal)
                if not data:
                    continue
                else:
                    fly_id += 1
            else:
                data = main(path, sex=sex, PARAMETERS = PARAMETERS, fly_id = fly_id, saccade=saccade, normal=normal)
                if not data:
                    continue
                else:
                    fly_id += 1
            if len(data) > 2:
                total_response_data = total_response_data + list(data[0])
                total_response_speed_data = total_response_speed_data + list(data[1])
                total_response_angular_speed_data = total_response_angular_speed_data + list(data[2])
                saccade_peak_list = saccade_peak_list + list(data[3][0])
                saccade_turn_list = saccade_turn_list + list(data[3][1])
                total_response_fwd_vel = total_response_fwd_vel + list(data[4])
                total_response_side_vel = total_response_side_vel + list(data[5])
                pos_list = pos_list + list(data[6])
                total_orientation_data = total_orientation_data + list(data[7])
                fly_name_dict['filename'].append(os.path.basename(data[8]['filename']))
                fly_name_dict['flyid'].append(data[8]['flyid'])
            else:
                if len(data) == 1:
                    continue
                else:
                    saccade_peak_list = saccade_peak_list + list(data[0])
                    saccade_turn_list = saccade_turn_list + list(data[1])
    # json.dump(fly_name_dict, flyname_file)
    fly_name_df = pd.DataFrame(fly_name_dict)
    fly_name_df.to_csv(flyname_file)
    if saccade=='N':
        pass
    else:
        saccade_peak_list = np.array(saccade_peak_list)
        saccade_turn_list = np.array(saccade_turn_list)
        output = r'{}/{}'.format(Dir, filename + '_opto_response_' + sex)
        np.save(output + '_saccade_peak', saccade_peak_list)
        np.save(output + '_saccade_turn', saccade_turn_list)
    if normal=='Y':
        total_response_data= np.array(total_response_data)
        total_orientation_data = np.array(total_orientation_data)
        total_response_speed_data = np.array(total_response_speed_data)
        total_response_angular_speed_data = np.array(total_response_angular_speed_data)
        total_response_fwd_vel = np.array(total_response_fwd_vel)
        total_response_side_vel = np.array(total_response_side_vel)
        pos_list = np.array(pos_list)

        output = r'{}/{}'.format(Dir, filename+'_opto_response_'+sex)
        np.save(output, total_response_data)
        np.save(output+'_speed', total_response_speed_data)
        np.save(output+'_angular_speed', total_response_angular_speed_data)
        np.save(output+'_fwd_vel', total_response_fwd_vel)
        np.save(output+'_side_vel', total_response_side_vel)
        np.save(output + '_pos', pos_list)
        np.save(output + '_orientation', total_orientation_data)
    else:
        pass
    return 1


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument("Directory", help="Directory to be processed. It should not sub-directories containing data for individual flies.\
                        If the directory is not present in the current working directory, the full path should be given.")
    parser.add_argument("-n", "--name", help="Genotype of fly", default=None, required=True, type=str)
    parser.add_argument("-p", "--parts", help="How many sections does the pinwheel stimulus have?", default=2, type=int, required=True)
    parser.add_argument("-c", "--combine", help="Whether to go through all the sub-directories, preprocess and combine the data. \
                        Search for .npy files. If found, do not do combine", default=False, action='store_true', dest='combine')
    parser.add_argument("-a", "--analysis", help="Whether to perform simple analysis. \
                        Explained more in readme.", default=True, action='store_true', dest='analysis')
    parser.add_argument("-s", "--saccade", help="Whether to perform saccade analysis. \
                        Explained more in readme.", default=True, action='store_true', dest='saccade')
    parser.add_argument("-no-a", "--no-analysis", help="Whether to perform simple analysis. \
                        Explained more in readme.", default=True, action='store_false', dest='analysis')
    parser.add_argument("-no-s", "--no-saccade", help="Whether to perform saccade analysis. \
                        Explained more in readme.", default=True, action='store_false', dest='saccade')

    args = parser.parse_args()
    print(args.combine, args.analysis, args.saccade)
    if os.path.exists(args.Directory):
        if not os.path.isdir(args.Directory):
            raise Exception('Please input a directory/folder, not a file')
        else:
            pass
    else:
        raise Exception('{} directory does not exist. Please check the the directory/folder you have entered'.format(args.Directory))

    if args.parts == 2:
        parameters = ['stim_left', 'stim_right']
        condn_indices = [1, 0]
    elif args.parts == 3:
        parameters = ['stim_left', 'stim_right', 'stim_center']
        condn_indices = [2, 1, 0]
    elif args.parts == 4:
        parameters = ['stim_left', 'stim_right', 'stim_center', 'stim_rear']
        condn_indices = [3, 2, 1, 0]
    else:
        raise Exception('Invalid number of parts in the the stimulus. Please enter 2, 3 or 4./nOr change the code in \
                        process_directory_for_optomotor_analysis.py to include more parts')

    if args.combine:
        folder_preprocessed = False
        for files in os.listdir(args.Directory):
            if files.endswith('.npy'):
                folder_preprocessed = True
                combine = input('It seems this directory might have already been pre-processed. /nAre you sure you want to continue, data will be over-written? (Y/N)')
                if combine.casefold() == 'Y'.casefold():
                    combine_data(args.Directory, filename=args.name, sex='', PARAMETERS=parameters, saccade='Y', normal='Y')
                break
        if not folder_preprocessed:
            combine_data(args.Directory, filename=args.name, sex='', PARAMETERS=parameters, saccade='Y', normal='Y')
    else:
        folder_preprocessed = False
        for files in os.listdir(args.Directory):
            if files.endswith('.npy'):
                folder_preprocessed = True
                break
        if not folder_preprocessed:
            combine = input('It seems the contents of this folder have not been preprocessed. \nDo you want to analyze it now? (Y/N)')
            if combine.casefold() == 'Y'.casefold():
                combine_data(args.Directory, filename=args.name, sex='', PARAMETERS=parameters, saccade='Y', normal='Y')

    if args.analysis:
        analysis(args.Directory, filename=args.name, sex='', p=args.parts, condn_indices=condn_indices, ignore_values=[], variables=[])

    if args.saccade:
        analyse_saccade(args.Directory, filename=args.name, sex='', p=args.parts, condn_indices=condn_indices, ignore_values=[], variables=[])