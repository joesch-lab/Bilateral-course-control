import pandas as pd
import os
from statannotations.Annotator import Annotator
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import itertools
from scipy.stats import ttest_1samp
from statannotations.stats.StatTest import StatTest


global all_color_data
all_color_data = loadmat('colorList.mat')['all_data']

color_dict = json.load(open('color_dict.json'))
shades_ND = create_color_gradient(color_dict['FlpND'][0], color_dict['FlpND'][-1], 20)
shades_D = create_color_gradient(color_dict['FlpD'][0], color_dict['FlpD'][-1], 20)
shades_DB331 = create_color_gradient(color_dict['FlpNDxDB331'][0], color_dict['FlpNDxDB331'][-1], 20)
shades_ShakB2 = create_color_gradient(color_dict['ShakB2'][0], color_dict['ShakB2'][-1], 20)
shades_CantonS = create_color_gradient((0.0, 0.0, 0.0), (1, 1, 1), 20)
shades_LPTCKir21 = create_color_gradient((0.95, 0.6, 0.0), (0.35, 0.1, 0.0), 20)
genotype_colors = {'CantonS_18':shades_CantonS[-1], 'LPTC_Kir21':shades_LPTCKir21[-1], 'FlpND_18':shades_ND[-1], 'FlpD_18':shades_D[-1],
                   'Kir21_control':shades_CantonS[-10], 'LPTC_control':shades_CantonS[-15], 'H2split_shibire_activated_all':[27/255,158/255,119/255],
                    'LPTC_tshgal80_shibire_activated':[117/255,112/255,179/255], 'shibire_temperature_control':[0.8, 0.8, 0.8],
                   'Kir21_tshgal80_emptygal4':[0.8, 0.8, 0.8]}

folders_list = ['monocular_stimulus.csv', 'hemistimulus_files.csv', 'hemistimulus_files_25.csv', 'hemistimulus_files_hemi_eyecovered.csv',\
                'hemistimulus_files_quarter.csv', 'hemistimulus_files_quarter4.csv']
folders_list = ['monocular_stimulus.csv','hemistimulus_files.csv']
genotype_folders = {}
for folders in folders_list:
    filename_folder = os.path.join(r'D:\Roshan\Project\PythonCodes\Codes\Plotting', folders)
    filedata = pd.read_csv(filename_folder, index_col=0)
    genotype_folders = {**genotype_folders, **filedata.to_dict()['destination']}
print(genotype_folders['FlpND_18'])
# genotype_folders = {'CantonS':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\CantonS_18\Hemisphere\Contrast_100\ProcessedData',
#                 'LPTC_Kir21':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\LPTC_Kir21\Hemisphere\Contrast_100\ProcessedData',
#                 'FlpND':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\FlpND_18\Hemisphere\Contrast_100\ProcessedData',
#                 'FlpD':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\FlpD_18\Hemisphere\Contrast_100\ProcessedData',
#                 'Kir21_control':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Kir21_control\Hemisphere\New\Contrast_100\ProcessedData',
#                 'LPTC_control':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\LPTC_Gal4_control\Hemisphere\Contrast_100\ProcessedData',
#                 'H2_Kir21':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\New_H2_Kir2.1_18\Hemisphere\Contrast_100\ProcessedData',
#                 'H2_TNT':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\H2_TNT_18\H2_TNT_18_new_new\Hemisphere\Contrast_100\ProcessedData',
#                 'H2_shi_control':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\H2split_shibire\Hemisphere\22C_control\New_22C\ProcessedData',
#                 'H2_shi_activated':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\H2split_shibire\Hemisphere\32C_activated\ProcessedData',
#                 'TNT_control':r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\TNT_control_18\Hemisphere\Contrast_100\LPTC_negative\ProcessedData'}


genotypes=[]

## Fig1 
def plot_smooth_saccadic_contribution(genotypes, to_plot = [], labels=[]):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    ## plot contribution of saccadic and smooth responses to the total response
    for filename in genotypes:
        Dir = genotype_folders[filename]
        print(filename, Dir)
        for files in os.listdir(Dir):
            print(files)
            if 'saccadic_smooth_response' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]

        fig, ax = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
        ## stacked plot
        if len(to_plot) == 0:
            to_plot = ['total_angular_speed', 'all_saccade_angular_speed', 'all_saccade_removed_angular_speed']
            labels = ['Total Response', 'Saccadic Response', 'Smooth Response']
        stim_types = ['FrontToBack', 'BackToFront', 'Full Rotation']
        # stim_types = ['Full Rotation', 'BackToFront', 'FrontToBack']
        color = ['lightsteelblue', 'lightpink']
        # color = [(241/255,163/255,64/255), (153/255,142/255,195/255)]
        # 153, 142, 195 = purple

        n = len(data[to_plot[0]])-1
        for i in range(len(data[to_plot[0]])):
            saccadic = savgol_filter(data['all_saccade_angular_speed'][i][181:] - np.mean(data['all_saccade_angular_speed'][i][181:300]), 31, 2)
            smooth = savgol_filter(data['all_saccade_removed_angular_speed'][i][181:] - np.mean(data['all_saccade_removed_angular_speed'][i][181:300]), 31, 2)
            total = savgol_filter(data['total_angular_speed'][i][181:] - np.mean(data['total_angular_speed'][i][181:300]), 31, 2)
            if (np.mean(saccadic[120:]) > 0 and np.mean(smooth[120:]) > 0) or (np.mean(saccadic[120:]) < 0 and np.mean(smooth[120:]) < 0):
                if np.abs(np.mean(smooth[120:])) < np.abs(np.mean(saccadic[120:])):
                    ax[n-i].fill_between(np.arange(smooth.shape[0]), np.zeros(smooth.shape[0]), smooth, color=(153/255,142/255,195/255), alpha=0.6, lw=0, label='Smooth response')
                    ax[n-i].fill_between(np.arange(smooth.shape[0]), smooth, total, color=(241/255,163/255,64/255), alpha=0.6, lw=0, label='Saccadic response')
                    ax[n-i].plot(total, lw=1.5, c='k', label=labels[0], alpha=0.8, solid_capstyle='butt')
                else:
                    ax[n-i].fill_between(np.arange(smooth.shape[0]), np.zeros(smooth.shape[0]), saccadic, color=(241/255,163/255,64/255), alpha=0.6, lw=0, label='Saccadic response')
                    ax[n-i].fill_between(np.arange(smooth.shape[0]), saccadic, total, color=(153/255,142/255,195/255), alpha=0.6, lw=0, label='Smooth response')
                    ax[n-i].plot(total, lw=1.5, c='k', label=labels[0], alpha=0.8, solid_capstyle='butt')
            else:
                if np.abs(np.mean(smooth[120:])) < np.abs(np.mean(saccadic[120:])):
                    ax[n - i].fill_between(np.arange(smooth.shape[0]), np.zeros(smooth.shape[0]), smooth, color=(153/255,142/255,195/255), alpha=0.8, lw=0, label='Smooth response')
                    ax[n - i].fill_between(np.arange(smooth.shape[0]), smooth, total, color=(241/255,163/255,64/255), alpha=0.6, lw=0, label='Saccadic response')
                    ax[n-i].plot(total, lw=1.5, c='k', label=labels[0], alpha=0.8, solid_capstyle='butt')
                else:
                    ax[n-i].fill_between(np.arange(smooth.shape[0]), np.zeros(smooth.shape[0]), saccadic, color=(241/255,163/255,64/255), alpha=0.6, lw=0, label='Saccadic response')
                    ax[n - i].fill_between(np.arange(smooth.shape[0]), saccadic, total, color=(153/255,142/255,195/255), alpha=0.68, lw=0, label='Smooth response')
                    ax[n-i].plot(total, lw=1.5, c='k', label=labels[0], alpha=0.8, solid_capstyle='butt')


        # n = len(data[to_plot[0]])-1
        # for i in range(len(data[to_plot[0]])):
        #     ## create an array for stackplot
        #     data_to_plot = np.zeros((len(to_plot), (2*stim_period)-180))
        #     for j in range(len(to_plot)-1):
        #         data_to_plot[j, :] = savgol_filter(data[to_plot[1+j]][i][181:] - np.mean(data[to_plot[1+j]][i][181:300]), 31,2)
        #     ax[n-i].plot(savgol_filter(data[to_plot[0]][i][181:] - np.mean(data[to_plot[0]][i][181:300]), 31, 2), lw=1.5, c='k', label=labels[0], alpha=0.8, solid_capstyle='butt')
        #     ax[n-i].stackplot(np.arange(0, data_to_plot.shape[1], 1),  data_to_plot, baseline='zero', labels=labels, colors=color, alpha=0.6)
            ax[n-i].axhline(0, 0.05, 0.95, ls=(0, (5, 10)), c='grey', lw=0.8, alpha=0.6)
            ax[n-i].axvline(x=120, ls=(0, (5, 10)), c='grey', lw=0.8, alpha=0.6)
            ax[n-i].set_xticks([0, 120, data[to_plot[0]][i].shape[0]-180])
            ax[n-i].set_xticklabels([str(round(x / 60)) for x in [0, 120, data[to_plot[0]][i].shape[0]-180]])
            ax[n-i].spines['bottom'].set_bounds(high=data[to_plot[0]][i].shape[0]-180, low=0)
            ax[n-i].set_title(stim_types[i])

        # ax[0].set_xticks([])
        # ax[1].set_xticks([])
        ax[0].set_xticks([0, 120, data[to_plot[0]][0].shape[0]-180])
        ax[0].set_xticklabels([str(round(x / 60)) for x in [0, 120, data[to_plot[0]][0].shape[0]-180]])
        ax[0].spines['bottom'].set_bounds(high=data[to_plot[0]][i].shape[0]-180, low=0)
        # ax[0].spines['bottom'].set_visible(False)
        # ax[1].spines['bottom'].set_visible(False)
        top = 150
        bottom = -100
        ax[0].set_ylim(top=top+top/10, bottom=bottom)
        ax[0].set_yticks([bottom, 0, top])
        ax[0].spines['left'].set_bounds(high=top, low=bottom)
        ax[1].spines['left'].set_bounds(high=top, low=bottom)
        ax[2].spines['left'].set_bounds(high=top, low=bottom)
        ax[2].legend()
        fig.savefig(os.path.join(Dir_output, filename+'saccades_comparison_turning.png'), dpi=300)
        fig.savefig(os.path.join(Dir_output, filename+'saccades_comparison_turning.pdf'), dpi=600, transparent=True)
        fig.savefig(os.path.join(Dir_output, filename + 'saccades_comparison_turning.svg'), dpi=600, transparent=True)
    return 1


def plot_saccade_numbers_and_strength(genotypes):
    # # saccade numbers (bar plots showing total number of saccades)
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]

    df = pd.DataFrame({'saccades':[0], 'saccade_type':['syn'], 'stim_type':['pre'], 'genotype':['control']})
    time = ''
    for Dir, current_genotype in zip(genotype_Dir, genotypes):
        for files in os.listdir(Dir):
            if 'timewise_saccade_hemi' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]
        labels = data['variables']

        ylim = [8, 0]
        data_to_plot_syn = []
        data_to_plot_anti = []

        for i in range(3):
            data_to_plot_syn += [data['syn_pre_saccade_dist'][i], data['syn_saccade_dist'][i]]
            data_to_plot_anti += [data['anti_pre_saccade_dist'][i], data['anti_saccade_dist'][i]]
            # data_to_plot_syn += [data['syn_saccade_dist'][i]]
            # data_to_plot_anti += [data['anti_saccade_dist'][i]]
            ## make dataframe
            df1 = pd.DataFrame({'saccades':data['syn_pre_saccade_dist'+time][i], 'saccade_type':['syn']*len(data['syn_pre_saccade_dist'+time][i]),
                                'stim_type':['pre']*len(data['syn_pre_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_pre_saccade_dist'+time][i])})
            df = pd.concat([df, df1], ignore_index=True)
            df2 = pd.DataFrame({'saccades':data['syn_saccade_dist'+time][i], 'saccade_type':['syn']*len(data['syn_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['syn_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_saccade_dist'+time][i])})
            df = pd.concat([df, df2], ignore_index=True)
            df3 = pd.DataFrame({'saccades':data['anti_pre_saccade_dist'+time][i], 'saccade_type':['anti']*len(data['anti_pre_saccade_dist'+time][i]),
                                'stim_type':['pre']*len(data['anti_pre_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_pre_saccade_dist'+time][i])})
            df = pd.concat([df, df3], ignore_index=True)
            df4 = pd.DataFrame({'saccades':data['anti_saccade_dist'+time][i], 'saccade_type':['anti']*len(data['anti_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_saccade_dist'+time][i])})
            df = pd.concat([df, df4], ignore_index=True)
            df5 = pd.DataFrame({'saccades':np.subtract(data['anti_saccade_dist'+time][i], data['anti_pre_saccade_dist'+time][i]), 'saccade_type':['anti_change']*len(data['anti_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_saccade_dist'+time][i])})
            df = pd.concat([df, df5], ignore_index=True)
            df6 = pd.DataFrame({'saccades':np.subtract(data['syn_saccade_dist'+time][i], data['syn_pre_saccade_dist'+time][i]), 'saccade_type':['syn_change']*len(data['syn_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['syn_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_saccade_dist'+time][i])})
            df = pd.concat([df, df6], ignore_index=True)
            df7 = pd.DataFrame({'saccades': np.abs(data['anti_saccade_turns'+time][i]), 'saccade_type': ['anti_turns']*len(data['anti_saccade_turns'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_turns'+time][i]), 'genotype': [current_genotype]*len(data['anti_saccade_turns'+time][i])})
            df = pd.concat([df, df7], ignore_index=True)
            df8 = pd.DataFrame({'saccades': data['syn_saccade_turns' + time][i], 'saccade_type': ['syn_turns'] * len(data['syn_saccade_turns' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_turns' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_turns' + time][i])})
            df = pd.concat([df, df8], ignore_index=True)
            df9 = pd.DataFrame({'saccades': data['syn_saccade_peaks' + time][i], 'saccade_type': ['syn_peaks'] * len(data['syn_saccade_peaks' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_peaks' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_peaks' + time][i])})
            df = pd.concat([df, df9], ignore_index=True)
            df10 = pd.DataFrame({'saccades': data['syn_saccade_peaks' + time][i], 'saccade_type': ['syn_peaks'] * len(data['syn_saccade_peaks' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_peaks' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_peaks' + time][i])})
            df = pd.concat([df, df10], ignore_index=True)

    cmap = make_color_pastel(ListedColormap(all_color_data['PiYG5']))
    color = [cmap(0.95), cmap(0.05)]
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    linewidths = [0.8]*len(genotypes)
    for i, current_genotype in enumerate(genotypes):
        sns.violinplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back', 'pre'], hue='saccade_type',
                       hue_order=['syn', 'anti'], split=True, palette=color, linewidth=linewidths[i], scale='width', cut=0, inner=None, ax = ax)

        sns.violinplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back', 'pre'], hue='saccade_type',
                       hue_order=['syn_turns', 'anti_turns'], split=True, palette=color, linewidth=linewidths[i], scale='width', cut=0, inner=None, ax = ax2)

        sns.violinplot(data=df[(df.genotype==current_genotype) & (df.stim_type!='pre')], y='saccades', x='saccade_type', order=['syn_turns', 'anti_turns'], hue='saccade_type',
                       hue_order=['syn_turns', 'anti_turns'], split=True, palette=color, linewidth=linewidths[i], scale='width', cut=0, inner=None, ax = ax3)

    ax.get_legend().remove()
    alphas = np.linspace(0.8, 0.4, len(genotypes))
    print(alphas)
    edge_color = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0, 0) for filename in genotypes]
    # for i, violin in enumerate(ax.collections):
    #     violin.set_alpha(alphas[i//6])
    #     violin.set_edgecolor(edge_color[i // 6])

    max=20
    min=0
    ax.set_ylim(top=max, bottom=min)
    ax.set_yticks(np.arange(min, max+1, 5))
    ax.spines['bottom'].set_position(('outward', 1))
    ax.set_ylabel('Number of saccades')
    ax.set_title(time)

    ax2.get_legend().remove()
    ax3.get_legend().remove()
    alphas = np.linspace(0.8, 0.4, len(genotypes))
    edge_color = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0, 0) for filename in genotypes]
    # for i, violin in enumerate(ax2.collections):
    #     violin.set_alpha(alphas[i//6])
    #     violin.set_edgecolor(edge_color[i // 6])
    ax2.set_ylim(top=90, bottom=0)
    ax2.spines['bottom'].set_position(('outward', 1))
    ax2.set_ylabel('Saccade Strength')
    ax2.set_title(time)

    ax3.set_ylim(top=90, bottom=0)
    ax3.spines['bottom'].set_position(('outward', 1))
    ax3.set_ylabel('Saccade Strength')
    ax3.set_title(time)

    filename = ''
    for gene in genotypes:
        filename = filename + gene
    fig.savefig(os.path.join(Dir_output, 'saccade_number_'+filename+time+'.png'))
    fig.savefig(os.path.join(Dir_output, 'saccade_number_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig.savefig(os.path.join(Dir_output, 'saccade_number_' + filename + time + '.svg'), transparent=True, dpi=600)

    fig2.savefig(os.path.join(Dir_output, 'saccade_strength_'+filename+time+'.png'))
    fig2.savefig(os.path.join(Dir_output, 'saccade_strength_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig2.savefig(os.path.join(Dir_output, 'saccade_strength_' + filename + time + '.svg'), transparent=True, dpi=600)

    fig3.savefig(os.path.join(Dir_output, 'saccade_strength_alltogether_'+filename+time+'.png'))
    fig3.savefig(os.path.join(Dir_output, 'saccade_strength_alltogether_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig3.savefig(os.path.join(Dir_output, 'saccade_strength_alltogether_' + filename + time + '.svg'), transparent=True, dpi=600)

    fig, ax = plt.subplots()
    linewidths = [0.8]*len(genotypes)
    for i, current_genotype in enumerate(genotypes):
        sns.violinplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Front to Back', 'Back to Front', 'Full rotation'], hue='saccade_type',
                       hue_order=['syn_change', 'anti_change'],
                       split=True, palette=color, linewidth=linewidths[i], scale='count', cut=0, inner=None)
    ax.get_legend().remove()
    alphas = np.linspace(0.8, 0.4, len(genotypes))
    edge_color = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0, 0) for filename in genotypes]
    for i, violin in enumerate(ax.collections):
        violin.set_alpha(alphas[i//6])
        violin.set_edgecolor(edge_color[i // 6])

    max=20
    min=-10
    ax.set_ylim(top=max, bottom=min)
    ax.set_yticks(np.arange(min, max+1, 5))
    # ax.set_title(filename)
    ax.axhline(y=0, lw=1, ls='--', c='grey')
    ax.set_ylabel('Change in saccades')
    ax.set_title(time)
    fig.savefig(os.path.join(Dir_output, 'saccade_number_change_'+filename+time+'.png'))
    fig.savefig(os.path.join(Dir_output, 'saccade_number_change_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig.savefig(os.path.join(Dir_output, 'saccade_number_change_' + filename + time + '.svg'), transparent=True, dpi=600)
    return 1


def plot_antisaccade_numbers_and_strength_per_fly(genotypes):
    # # saccade numbers (bar plots showing total number of saccades)
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures\NewPlotsSuppl2024'
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]

    df = pd.DataFrame({'saccades':[0], 'saccade_type':['syn'], 'stim_type':['pre'], 'genotype':['control']})
    time = ''
    for Dir, current_genotype in zip(genotype_Dir, genotypes):
        for files in os.listdir(Dir):
            if 'timewise_saccade_hemi_per_fly' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]
        labels = data['variables']

        for i in range(3):
            ## make dataframe
            df1 = pd.DataFrame({'saccades':data['syn_pre_saccade_dist'+time][i], 'saccade_type':['syn']*len(data['syn_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['syn_pre_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_pre_saccade_dist'+time][i])})
            df = pd.concat([df, df1], ignore_index=True)
            df2 = pd.DataFrame({'saccades':data['syn_saccade_dist'+time][i], 'saccade_type':['syn']*len(data['syn_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['syn_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_saccade_dist'+time][i])})
            df = pd.concat([df, df2], ignore_index=True)
            df3 = pd.DataFrame({'saccades':data['anti_pre_saccade_dist'+time][i], 'saccade_type':['anti']*len(data['anti_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['anti_pre_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_pre_saccade_dist'+time][i])})
            df = pd.concat([df, df3], ignore_index=True)
            df4 = pd.DataFrame({'saccades':data['anti_saccade_dist'+time][i], 'saccade_type':['anti']*len(data['anti_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_saccade_dist'+time][i])})
            df = pd.concat([df, df4], ignore_index=True)
            df5 = pd.DataFrame({'saccades':np.subtract(data['anti_saccade_dist'+time][i], data['anti_pre_saccade_dist'+time][i]), 'saccade_type':['anti_change']*len(data['anti_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_saccade_dist'+time][i])})
            df = pd.concat([df, df5], ignore_index=True)
            df6 = pd.DataFrame({'saccades':np.subtract(data['syn_saccade_dist'+time][i], data['syn_pre_saccade_dist'+time][i]), 'saccade_type':['syn_change']*len(data['syn_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['syn_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_saccade_dist'+time][i])})
            df = pd.concat([df, df6], ignore_index=True)
            df7 = pd.DataFrame({'saccades': np.abs(data['anti_saccade_turns'+time][i]), 'saccade_type': ['anti_turns']*len(data['anti_saccade_turns'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_turns'+time][i]), 'genotype': [current_genotype]*len(data['anti_saccade_turns'+time][i])})
            df = pd.concat([df, df7], ignore_index=True)
            df8 = pd.DataFrame({'saccades': data['syn_saccade_turns' + time][i], 'saccade_type': ['syn_turns'] * len(data['syn_saccade_turns' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_turns' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_turns' + time][i])})
            df = pd.concat([df, df8], ignore_index=True)

            df9 = pd.DataFrame({'saccades': np.abs(data['anti_saccade_peaks'+time][i]), 'saccade_type': ['anti_peaks']*len(data['anti_saccade_peaks'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_peaks'+time][i]), 'genotype': [current_genotype]*len(data['anti_saccade_peaks'+time][i])})
            df = pd.concat([df, df9], ignore_index=True)
            df10 = pd.DataFrame({'saccades': data['syn_saccade_peaks' + time][i], 'saccade_type': ['syn_peaks'] * len(data['syn_saccade_peaks' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_peaks' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_peaks' + time][i])})
            df = pd.concat([df, df10], ignore_index=True)

            df11 = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_turns'+time][i]), 'saccade_type': ['pre_anti_turns']*len(data['pre_anti_saccade_turns'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['pre_anti_saccade_turns'+time][i]), 'genotype': [current_genotype]*len(data['pre_anti_saccade_turns'+time][i])})
            df = pd.concat([df, df11], ignore_index=True)
            df12 = pd.DataFrame({'saccades': data['pre_syn_saccade_turns' + time][i], 'saccade_type': ['pre_syn_turns'] * len(data['pre_syn_saccade_turns' + time][i]),
                                'stim_type': [labels[i]+'pre'] * len(data['pre_syn_saccade_turns' + time][i]), 'genotype': [current_genotype] * len(data['pre_syn_saccade_turns' + time][i])})
            df = pd.concat([df, df12], ignore_index=True)

            df12 = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_peaks'+time][i]), 'saccade_type': ['pre_anti_peaks']*len(data['pre_anti_saccade_peaks'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['pre_anti_saccade_peaks'+time][i]), 'genotype': [current_genotype]*len(data['pre_anti_saccade_peaks'+time][i])})
            df = pd.concat([df, df12], ignore_index=True)
            df13 = pd.DataFrame({'saccades': data['pre_syn_saccade_peaks' + time][i], 'saccade_type': ['pre_syn_peaks'] * len(data['pre_syn_saccade_peaks' + time][i]),
                                'stim_type': [labels[i]+'pre'] * len(data['pre_syn_saccade_peaks' + time][i]), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks' + time][i])})
            df = pd.concat([df, df13], ignore_index=True)

        df_new = pd.DataFrame({'saccades': data['syn_saccade_turns_allstimuli'], 'saccade_type': ['syn_turns'] * len(data['syn_saccade_turns_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['syn_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['syn_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': data['pre_syn_saccade_turns_allstimuli'], 'saccade_type': ['pre_syn_turns'] * len(data['pre_syn_saccade_turns_allstimuli']),
                        'stim_type': ['all_pre'] * len(data['pre_syn_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': data['syn_saccade_peaks_allstimuli'], 'saccade_type': ['syn_peaks'] * len(data['syn_saccade_peaks_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': data['pre_syn_saccade_peaks_allstimuli'], 'saccade_type': ['pre_syn_peaks'] * len(data['pre_syn_saccade_peaks_allstimuli']),
                        'stim_type': ['all_pre'] * len(data['pre_syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.nanmean([data['pre_syn_saccade_peaks_allstimuli'], data['syn_saccade_peaks_allstimuli']], axis=0), 'saccade_type': ['all_syn_peaks'] * len(data['pre_syn_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.nanmean([data['pre_syn_saccade_turns_allstimuli'], data['syn_saccade_turns_allstimuli']], axis=0), 'saccade_type': ['all_syn_turns'] * len(data['pre_syn_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)

        df_new = pd.DataFrame({'saccades': np.abs(data['anti_saccade_turns_allstimuli']), 'saccade_type': ['anti_turns'] * len(data['anti_saccade_turns_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['anti_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['anti_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_turns_allstimuli']), 'saccade_type': ['pre_anti_turns'] * len(data['pre_anti_saccade_turns_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['pre_anti_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(data['anti_saccade_peaks_allstimuli']), 'saccade_type': ['anti_peaks'] * len(data['anti_saccade_peaks_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_peaks_allstimuli']), 'saccade_type': ['anti_syn_peaks'] * len(data['pre_anti_saccade_peaks_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['pre_anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(np.nanmean([data['pre_anti_saccade_peaks_allstimuli'], data['anti_saccade_peaks_allstimuli']], axis=0)), 'saccade_type': ['all_anti_peaks'] * len(data['pre_anti_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(np.nanmean([data['pre_anti_saccade_turns_allstimuli'], data['anti_saccade_turns_allstimuli']], axis=0)), 'saccade_type': ['all_anti_turns'] * len(data['pre_anti_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)

    names = ''
    for gen in genotypes:
        names+=gen

    cmap = make_color_pastel(ListedColormap(all_color_data['PiYG5']))
    color = [cmap(0.95), cmap(0.05)]

    # fig, ax = plt.subplots(figsize=(1*len(genotypes), 3+1))
    # sns.pointplot(data=df[df.stim_type == 'pre_and_stim'][df.saccade_type == 'all_anti_turns'], y='saccades', x='genotype',
    #               order=genotypes, dodge=0.4, ax=ax, palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
    # plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
    # plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
    # sns.stripplot(data=df[df.stim_type == 'pre_and_stim'][df.saccade_type == 'all_anti_turns'], y='saccades', x='genotype',
    #               order=genotypes, palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax=ax)
    #
    # pairs = list(itertools.combinations(list(itertools.product(genotypes, labels, repeat=1)), 2))
    # print(pairs)
    # annotator = Annotator(data=df[df.stim_type == 'pre_and_stim'][df.saccade_type == 'all_anti_turns'], y='saccades', x='genotype',
    #                       order=genotypes, pairs=pairs, ax=ax)
    # annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    # _, annotations = annotator.apply_and_annotate()
    # save_stats_data(annotations, Dir_output, names + '_antisaccade_strength_allstimandpre')
    # ax.get_legend().remove()
    # ax.set_ylim(top=65, bottom=0)
    # ax.spines['bottom'].set_position(('outward', 1))
    # ax.set_ylabel('Saccade Strength')
    # ax.set_title(time)
    #
    # fig.savefig(os.path.join(Dir_output, names+'antisaccade_strength_allstimandpre.png'))
    # fig.savefig(os.path.join(Dir_output, names+'antisaccade_strength_allstimandpre.pdf'), dpi=600, transparent=True)
    #
    # fig, ax = plt.subplots(figsize=(1*len(genotypes), 3+1))
    # sns.pointplot(data=df[df.stim_type == 'pre_and_stim'][df.saccade_type == 'all_anti_peaks'], y='saccades', x='genotype', order=genotypes,
    #               dodge=0.4, ax=ax, palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
    # plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
    # plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
    # sns.stripplot(data=df[df.stim_type == 'pre_and_stim'][df.saccade_type == 'all_anti_peaks'], y='saccades', x='genotype',
    #               order=genotypes, palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax=ax)
    #
    # pairs = list(itertools.combinations(list(itertools.product(genotypes,labels, repeat=1)), 2))
    # annotator = Annotator(data=df[df.stim_type == 'pre_and_stim'][df.saccade_type == 'all_anti_peaks'], y='saccades', x='genotype',
    #                       order=genotypes, pairs=pairs, ax=ax)
    # annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    # _, annotations = annotator.apply_and_annotate()
    # save_stats_data(annotations, Dir_output, names + '_antisaccade_peak_allstimandpre')
    #
    # ax.get_legend().remove()
    # ax.set_ylim(top=750, bottom=250)
    # ax.spines['bottom'].set_position(('outward', 1))
    # ax.set_ylabel('Saccade Strength')
    # ax.set_title(time)
    #
    # fig.savefig(os.path.join(Dir_output, names+'antisaccade_peak_allstimandpre.png'))
    # fig.savefig(os.path.join(Dir_output, names+'antisaccade_peak_allstimandpre.pdf'), dpi=600, transparent=True)

    fig, ax = plt.subplots(1, len(genotypes), figsize=(3*len(genotypes), 3))
    fig2, ax2 = plt.subplots(1, len(genotypes), figsize=(3*len(genotypes), 3))
    palette = [(0, 0, 0) for filename in genotypes]
    for j, current_genotype in enumerate(genotypes):
        if len(genotypes) == 1:
            axis=ax
            axis2=ax2
        else:
            axis=ax[j]
            axis2=ax2[j]
        sns.pointplot(data=df[df.genotype==current_genotype][df.saccade_type=='anti'], y='saccades', x='stim_type',
                    order=['Full rotation', 'Back to Front', 'Front to Back'], dodge=0.4, ax=axis, palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(axis.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(axis.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(data=df[df.genotype==current_genotype][df.saccade_type=='anti'], y='saccades', x='stim_type',
                      order=['Full rotation', 'Back to Front', 'Front to Back'], palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=False, ax = axis)

        axis.set_title(current_genotype)
        pairs = list(itertools.combinations(['Full rotation', 'Back to Front', 'Front to Back'], 2))
        annotator = Annotator(axis, data=df[df.genotype==current_genotype][df.saccade_type=='anti'], y='saccades', x='stim_type',
                              order=['Full rotation', 'Back to Front', 'Front to Back'], pairs=pairs)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, current_genotype + '_antisaccade_number_per_fly')

        max = 15
        min = 0
        axis.set_ylim(top=max, bottom=min-2)
        axis.set_yticks(np.arange(min, max+1, 5))
        axis.spines['bottom'].set_position(('outward', 1))
        axis.spines['left'].set_bounds(low=min, high=max)
        axis.set_ylabel('Number of saccades')
        axis.set_title(time)

        sns.pointplot(data=df[df.genotype==current_genotype][df.saccade_type=='anti_turns'], y='saccades', x='stim_type',
                      order=['Full rotation', 'Back to Front', 'Front to Back'], dodge=0.4, ax=axis2, palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(axis2.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(axis2.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(data=df[df.genotype==current_genotype][df.saccade_type=='anti_turns'], y='saccades', x='stim_type',
                      order=['Full rotation', 'Back to Front', 'Front to Back'], palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=False, ax=axis2)

        axis2.set_title(current_genotype)
        pairs = list(itertools.combinations(['Full rotation', 'Back to Front', 'Front to Back'], 2))
        annotator = Annotator(axis2, data=df[df.genotype==current_genotype][df.saccade_type=='anti_turns'], y='saccades', x='stim_type',
                              order=['Full rotation', 'Back to Front', 'Front to Back'], pairs=pairs)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, current_genotype + '_antisaccade_strength_per_fly')

        # axis2.get_legend().remove()
        axis2.set_ylim(top=75, bottom=0)
        # ax2.set_yticks(np.arange(min, max+1, 5))
        axis2.spines['bottom'].set_position(('outward', 1))
        axis2.set_ylabel('Saccade Strength')
        axis2.set_title(time)

    filename = ''
    for gene in genotypes:
        filename = filename + gene

    fig.savefig(os.path.join(Dir_output, 'antisaccade_number_per_fly_'+filename+time+'.png'))
    fig.savefig(os.path.join(Dir_output, 'antisaccade_number_per_fly_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig.savefig(os.path.join(Dir_output, 'antisaccade_number_per_fly_' + filename + time + '.svg'), transparent=True, dpi=600)

    fig2.savefig(os.path.join(Dir_output, 'antisaccade_strength_per_fly_'+filename+time+'.png'))
    fig2.savefig(os.path.join(Dir_output, 'antisaccade_strength_per_fly_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig2.savefig(os.path.join(Dir_output, 'antisaccade_strength_per_fly_' + filename + time + '.svg'), transparent=True, dpi=600)

    # fig, ax = plt.subplots(1, len(genotypes), figsize=(3*len(genotypes), 3))
    # print(len(genotypes))
    # for j, current_genotype in enumerate(genotypes):
    #     sns.pointplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
    #                 hue_order=['syn_change', 'anti_change'], dodge=0.4, ax=ax[j], palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
    #     plt.setp(ax[j].lines, zorder=100, alpha=0.95, lw=1.5)
    #     plt.setp(ax[j].collections, zorder=100, alpha=0.95, lw=1.5)
    #     sns.stripplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
    #                    hue_order=['syn_change', 'anti_change'], palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax = ax[j])
    #
    #     ax[j].set_title(current_genotype)
    #     pairs = []
    #     for x in ['Full rotation', 'Back to Front', 'Front to Back']:
    #         pairs.append(((x, 'syn_change'), (x, 'anti_change')))
    #     annotator = Annotator(ax[j], data=df[df.genotype==current_genotype], y='saccades', x='stim_type', hue='saccade_type', hue_order=['syn_change', 'anti_change'],
    #                           order=['Full rotation', 'Back to Front', 'Front to Back'], pairs=pairs)
    #     annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    #     _, annotations = annotator.apply_and_annotate()
    #     save_stats_data(annotations, Dir_output, names + '_saccade_number_change_per_fly')
    #
    # max=10
    # min=-5
    # for i in range(len(genotypes)):
    #     ax[i].get_legend().remove()
    #     ax[i].set_ylim(top=max, bottom=min)
    #     ax[i].set_yticks(np.arange(min, max+1, 5))
    #     # ax.set_title(filename)
    #     ax[i].axhline(y=0, lw=1, ls='--', c='grey')
    #     ax[i].set_ylabel('Change in saccades')
    #     ax[i].set_title(time)
    # fig.savefig(os.path.join(Dir_output, 'saccade_number_change_per_fly_'+filename+time+'.png'))
    # fig.savefig(os.path.join(Dir_output, 'saccade_number_change_per_fly_'+filename+time+'.pdf'), transparent=True, dpi=600)
    # fig.savefig(os.path.join(Dir_output, 'saccade_number_change_per_fly_' + filename + time + '.svg'), transparent=True, dpi=600)
    return 1


def plot_saccade_numbers_and_strength_per_fly(genotypes):
    # # saccade numbers (bar plots showing total number of saccades)
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures\NewPlotsSuppl2024'
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]

    df = pd.DataFrame({'saccades':[0], 'saccade_type':['syn'], 'stim_type':['pre'], 'genotype':['control']})
    time = ''
    for Dir, current_genotype in zip(genotype_Dir, genotypes):
        for files in os.listdir(Dir):
            if 'timewise_saccade_hemi_per_fly' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]
        labels = data['variables']

        for i in range(3):
            ## make dataframe
            df1 = pd.DataFrame({'saccades':data['syn_pre_saccade_dist'+time][i], 'saccade_type':['syn']*len(data['syn_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['syn_pre_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_pre_saccade_dist'+time][i])})
            df = pd.concat([df, df1], ignore_index=True)
            df2 = pd.DataFrame({'saccades':data['syn_saccade_dist'+time][i], 'saccade_type':['syn']*len(data['syn_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['syn_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_saccade_dist'+time][i])})
            df = pd.concat([df, df2], ignore_index=True)
            df3 = pd.DataFrame({'saccades':data['anti_pre_saccade_dist'+time][i], 'saccade_type':['anti']*len(data['anti_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['anti_pre_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_pre_saccade_dist'+time][i])})
            df = pd.concat([df, df3], ignore_index=True)
            df4 = pd.DataFrame({'saccades':data['anti_saccade_dist'+time][i], 'saccade_type':['anti']*len(data['anti_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_saccade_dist'+time][i])})
            df = pd.concat([df, df4], ignore_index=True)
            df5 = pd.DataFrame({'saccades':np.subtract(data['anti_saccade_dist'+time][i], data['anti_pre_saccade_dist'+time][i]), 'saccade_type':['anti_change']*len(data['anti_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['anti_saccade_dist'+time][i])})
            df = pd.concat([df, df5], ignore_index=True)
            df6 = pd.DataFrame({'saccades':np.subtract(data['syn_saccade_dist'+time][i], data['syn_pre_saccade_dist'+time][i]), 'saccade_type':['syn_change']*len(data['syn_pre_saccade_dist'+time][i]),
                                'stim_type':[labels[i]]*len(data['syn_saccade_dist'+time][i]), 'genotype':[current_genotype]*len(data['syn_saccade_dist'+time][i])})
            df = pd.concat([df, df6], ignore_index=True)
            df7 = pd.DataFrame({'saccades': np.abs(data['anti_saccade_turns'+time][i]), 'saccade_type': ['anti_turns']*len(data['anti_saccade_turns'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_turns'+time][i]), 'genotype': [current_genotype]*len(data['anti_saccade_turns'+time][i])})
            df = pd.concat([df, df7], ignore_index=True)
            df8 = pd.DataFrame({'saccades': data['syn_saccade_turns' + time][i], 'saccade_type': ['syn_turns'] * len(data['syn_saccade_turns' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_turns' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_turns' + time][i])})
            df = pd.concat([df, df8], ignore_index=True)

            df9 = pd.DataFrame({'saccades': np.abs(data['anti_saccade_peaks'+time][i]), 'saccade_type': ['anti_peaks']*len(data['anti_saccade_peaks'+time][i]),
                                'stim_type':[labels[i]]*len(data['anti_saccade_peaks'+time][i]), 'genotype': [current_genotype]*len(data['anti_saccade_peaks'+time][i])})
            df = pd.concat([df, df9], ignore_index=True)
            df10 = pd.DataFrame({'saccades': data['syn_saccade_peaks' + time][i], 'saccade_type': ['syn_peaks'] * len(data['syn_saccade_peaks' + time][i]),
                                'stim_type': [labels[i]] * len(data['syn_saccade_peaks' + time][i]), 'genotype': [current_genotype] * len(data['syn_saccade_peaks' + time][i])})
            df = pd.concat([df, df10], ignore_index=True)

            df11 = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_turns'+time][i]), 'saccade_type': ['pre_anti_turns']*len(data['pre_anti_saccade_turns'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['pre_anti_saccade_turns'+time][i]), 'genotype': [current_genotype]*len(data['pre_anti_saccade_turns'+time][i])})
            df = pd.concat([df, df11], ignore_index=True)
            df12 = pd.DataFrame({'saccades': data['pre_syn_saccade_turns' + time][i], 'saccade_type': ['pre_syn_turns'] * len(data['pre_syn_saccade_turns' + time][i]),
                                'stim_type': [labels[i]+'pre'] * len(data['pre_syn_saccade_turns' + time][i]), 'genotype': [current_genotype] * len(data['pre_syn_saccade_turns' + time][i])})
            df = pd.concat([df, df12], ignore_index=True)

            df12 = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_peaks'+time][i]), 'saccade_type': ['pre_anti_peaks']*len(data['pre_anti_saccade_peaks'+time][i]),
                                'stim_type':[labels[i]+'pre']*len(data['pre_anti_saccade_peaks'+time][i]), 'genotype': [current_genotype]*len(data['pre_anti_saccade_peaks'+time][i])})
            df = pd.concat([df, df12], ignore_index=True)
            df13 = pd.DataFrame({'saccades': data['pre_syn_saccade_peaks' + time][i], 'saccade_type': ['pre_syn_peaks'] * len(data['pre_syn_saccade_peaks' + time][i]),
                                'stim_type': [labels[i]+'pre'] * len(data['pre_syn_saccade_peaks' + time][i]), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks' + time][i])})
            df = pd.concat([df, df13], ignore_index=True)

        df_new = pd.DataFrame({'saccades': data['syn_saccade_turns_allstimuli'], 'saccade_type': ['syn_turns'] * len(data['syn_saccade_turns_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['syn_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['syn_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': data['pre_syn_saccade_turns_allstimuli'], 'saccade_type': ['pre_syn_turns'] * len(data['pre_syn_saccade_turns_allstimuli']),
                        'stim_type': ['all_pre'] * len(data['pre_syn_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': data['syn_saccade_peaks_allstimuli'], 'saccade_type': ['syn_peaks'] * len(data['syn_saccade_peaks_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': data['pre_syn_saccade_peaks_allstimuli'], 'saccade_type': ['pre_syn_peaks'] * len(data['pre_syn_saccade_peaks_allstimuli']),
                        'stim_type': ['all_pre'] * len(data['pre_syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.nanmean([data['pre_syn_saccade_peaks_allstimuli'], data['syn_saccade_peaks_allstimuli']], axis=0), 'saccade_type': ['all_syn_peaks'] * len(data['pre_syn_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.nanmean([data['pre_syn_saccade_turns_allstimuli'], data['syn_saccade_turns_allstimuli']], axis=0), 'saccade_type': ['all_syn_turns'] * len(data['pre_syn_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_syn_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_syn_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)

        df_new = pd.DataFrame({'saccades': np.abs(data['anti_saccade_turns_allstimuli']), 'saccade_type': ['anti_turns'] * len(data['anti_saccade_turns_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['anti_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['anti_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_turns_allstimuli']), 'saccade_type': ['pre_anti_turns'] * len(data['pre_anti_saccade_turns_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['pre_anti_saccade_turns_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_turns_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(data['anti_saccade_peaks_allstimuli']), 'saccade_type': ['anti_peaks'] * len(data['anti_saccade_peaks_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(data['pre_anti_saccade_peaks_allstimuli']), 'saccade_type': ['anti_syn_peaks'] * len(data['pre_anti_saccade_peaks_allstimuli']),
                        'stim_type': ['all_stimuli'] * len(data['pre_anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(np.nanmean([data['pre_anti_saccade_peaks_allstimuli'], data['anti_saccade_peaks_allstimuli']], axis=0)), 'saccade_type': ['all_anti_peaks'] * len(data['pre_anti_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)
        df_new = pd.DataFrame({'saccades': np.abs(np.nanmean([data['pre_anti_saccade_turns_allstimuli'], data['anti_saccade_turns_allstimuli']], axis=0)), 'saccade_type': ['all_anti_turns'] * len(data['pre_anti_saccade_peaks_allstimuli']),
                        'stim_type': ['pre_and_stim'] * len(data['pre_anti_saccade_peaks_allstimuli']), 'genotype': [current_genotype] * len(data['pre_anti_saccade_peaks_allstimuli'])})
        df = pd.concat([df, df_new], ignore_index=True)

    names = ''
    for gen in genotypes:
        names+=gen
    # for saccade_type in ['syn', 'anti']:
    #     for stim_type in ['Full rotation', 'Back to Front', 'Front to Back']:
    #         print(stim_type, saccade_type)
    #         U, p = mannwhitneyu(df[(df.genotype==genotypes[0]) & (df.saccade_type==saccade_type) & (df.stim_type==stim_type)]['saccades'],
    #                             df[(df.genotype==genotypes[1]) & (df.saccade_type==saccade_type) & (df.stim_type==stim_type)]['saccades'], nan_policy='omit')
    #         print(U, p)
    # return 1

    cmap = make_color_pastel(ListedColormap(all_color_data['PiYG5']))
    color = [cmap(0.95), cmap(0.05)]

    fig, ax = plt.subplots(figsize=(1*len(genotypes), 3+1))
    sns.pointplot(data=df[df.stim_type == 'pre_and_stim'], y='saccades', x='genotype', order=genotypes, hue='saccade_type',
                  hue_order=['all_syn_turns', 'all_anti_turns'], dodge=0.4, ax=ax, palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
    plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
    plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
    sns.stripplot(data=df[df.stim_type == 'pre_and_stim'], y='saccades', x='genotype', order=genotypes, hue='saccade_type',
                  hue_order=['all_syn_turns', 'all_anti_turns'], palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax=ax)
    pairs = []
    for x in genotypes:
        pairs.append(((x, 'all_syn_turns'), (x, 'all_anti_turns')))
    for i in range(1,len(genotypes)):
        for y in ['all_syn_turns', 'all_anti_turns']:
            pairs.append(((genotypes[0], y), (genotypes[i], y)))
    annotator = Annotator(data=df[df.stim_type == 'pre_and_stim'], y='saccades', x='genotype', order=genotypes, hue='saccade_type',
                  hue_order=['all_syn_turns', 'all_anti_turns'], pairs=pairs, ax=ax)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    _, annotations = annotator.apply_and_annotate()
    save_stats_data(annotations, Dir_output, names + '_saccade_strength_allstimandpre')
    ax.get_legend().remove()
    ax.set_ylim(top=65, bottom=0)
    ax.spines['bottom'].set_position(('outward', 1))
    ax.set_ylabel('Saccade Strength')
    ax.set_title(time)

    fig.savefig(os.path.join(Dir_output, names+'saccade_strength_allstimandpre.png'))
    fig.savefig(os.path.join(Dir_output, names+'saccade_strength_allstimandpre.pdf'), dpi=600, transparent=True)

    fig, ax = plt.subplots(figsize=(1*len(genotypes), 3+1))
    sns.pointplot(data=df[df.stim_type == 'pre_and_stim'], y='saccades', x='genotype', order=genotypes, hue='saccade_type',
                  hue_order=['all_syn_peaks', 'all_anti_peaks'], dodge=0.4, ax=ax, palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
    plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
    plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
    sns.stripplot(data=df[df.stim_type == 'pre_and_stim'], y='saccades', x='genotype', order=genotypes, hue='saccade_type',
                  hue_order=['all_syn_peaks', 'all_anti_peaks'], palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax=ax)
    pairs = []
    for x in genotypes:
        pairs.append(((x, 'all_syn_peaks'), (x, 'all_anti_peaks')))
    for i in range(1,len(genotypes)):
        for y in ['all_syn_peaks', 'all_anti_peaks']:
            pairs.append(((genotypes[0], y), (genotypes[i], y)))
    annotator = Annotator(data=df[df.stim_type == 'pre_and_stim'], y='saccades', x='genotype', order=genotypes, hue='saccade_type',
                  hue_order=['all_syn_peaks', 'all_anti_peaks'], pairs=pairs, ax=ax)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    _, annotations = annotator.apply_and_annotate()
    save_stats_data(annotations, Dir_output, names + '_saccade_peak_allstimandpre')

    ax.get_legend().remove()
    ax.set_ylim(top=750, bottom=250)
    ax.spines['bottom'].set_position(('outward', 1))
    ax.set_ylabel('Saccade Strength')
    ax.set_title(time)

    fig.savefig(os.path.join(Dir_output, names+'saccade_peak_allstimandpre.png'))
    fig.savefig(os.path.join(Dir_output, names+'saccade_peak_allstimandpre.pdf'), dpi=600, transparent=True)

    fig, ax = plt.subplots(1, len(genotypes), figsize=(3*len(genotypes), 3))
    fig2, ax2 = plt.subplots(1, len(genotypes), figsize=(3*len(genotypes), 3))

    for j, current_genotype in enumerate(genotypes):
        if len(genotypes) == 1:
            axis=ax
            axis2=ax2
        else:
            axis=ax[j]
            axis2=ax2[j]
        sns.pointplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
                    hue_order=['syn', 'anti'], dodge=0.4, ax=axis, palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(axis.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(axis.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
                       hue_order=['syn', 'anti'], palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax = axis)

        axis.set_title(current_genotype)
        pairs = []
        for x in ['Full rotation', 'Back to Front', 'Front to Back']:
            pairs.append(((x, 'syn'), (x, 'anti')))
        annotator = Annotator(axis, data=df[df.genotype==current_genotype], y='saccades', x='stim_type', hue='saccade_type', hue_order=['syn', 'anti'], order=['Full rotation', 'Back to Front', 'Front to Back'], pairs=pairs)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, current_genotype + '_saccade_number_per_fly')

        max = 15
        min = 0
        axis.set_ylim(top=max, bottom=min-2)
        axis.set_yticks(np.arange(min, max+1, 5))
        axis.spines['bottom'].set_position(('outward', 1))
        axis.spines['left'].set_bounds(low=min, high=max)
        axis.set_ylabel('Number of saccades')
        axis.set_title(time)


        sns.pointplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
                    hue_order=['syn_turns', 'anti_turns'], dodge=0.4, ax=axis2, palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(axis2.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(axis2.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
                       hue_order=['syn_turns', 'anti_turns'], palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax=axis2)

        axis2.set_title(current_genotype)
        pairs = []
        for x in ['Full rotation', 'Back to Front', 'Front to Back']:
            pairs.append(((x, 'syn_turns'), (x, 'anti_turns')))
        annotator = Annotator(axis2, data=df[df.genotype==current_genotype], y='saccades', x='stim_type', hue='saccade_type', hue_order=['syn_turns', 'anti_turns'],
                              order=['Full rotation', 'Back to Front', 'Front to Back'], pairs=pairs)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, current_genotype + '_saccade_strength_per_fly')

        axis2.get_legend().remove()
        axis2.set_ylim(top=75, bottom=0)
        # ax2.set_yticks(np.arange(min, max+1, 5))
        axis2.spines['bottom'].set_position(('outward', 1))
        axis2.set_ylabel('Saccade Strength')
        axis2.set_title(time)

    filename = ''
    for gene in genotypes:
        filename = filename + gene

    fig.savefig(os.path.join(Dir_output, 'saccade_number_per_fly_'+filename+time+'.png'))
    fig.savefig(os.path.join(Dir_output, 'saccade_number_per_fly_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig.savefig(os.path.join(Dir_output, 'saccade_number_per_fly_' + filename + time + '.svg'), transparent=True, dpi=600)

    fig2.savefig(os.path.join(Dir_output, 'saccade_strength_per_fly_'+filename+time+'.png'))
    fig2.savefig(os.path.join(Dir_output, 'saccade_strength_per_fly_'+filename+time+'.pdf'), transparent=True, dpi=600)
    fig2.savefig(os.path.join(Dir_output, 'saccade_strength_per_fly_' + filename + time + '.svg'), transparent=True, dpi=600)

    # fig, ax = plt.subplots(1, len(genotypes), figsize=(3*len(genotypes), 3))
    # print(len(genotypes))
    # for j, current_genotype in enumerate(genotypes):
    #     sns.pointplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
    #                 hue_order=['syn_change', 'anti_change'], dodge=0.4, ax=ax[j], palette=color, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
    #     plt.setp(ax[j].lines, zorder=100, alpha=0.95, lw=1.5)
    #     plt.setp(ax[j].collections, zorder=100, alpha=0.95, lw=1.5)
    #     sns.stripplot(data=df[df.genotype==current_genotype], y='saccades', x='stim_type', order=['Full rotation', 'Back to Front', 'Front to Back'], hue='saccade_type',
    #                    hue_order=['syn_change', 'anti_change'], palette=color, size=7, jitter=0.15, linewidth=0.0, dodge=True, ax = ax[j])
    #
    #     ax[j].set_title(current_genotype)
    #     pairs = []
    #     for x in ['Full rotation', 'Back to Front', 'Front to Back']:
    #         pairs.append(((x, 'syn_change'), (x, 'anti_change')))
    #     annotator = Annotator(ax[j], data=df[df.genotype==current_genotype], y='saccades', x='stim_type', hue='saccade_type', hue_order=['syn_change', 'anti_change'],
    #                           order=['Full rotation', 'Back to Front', 'Front to Back'], pairs=pairs)
    #     annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    #     _, annotations = annotator.apply_and_annotate()
    #     save_stats_data(annotations, Dir_output, names + '_saccade_number_change_per_fly')
    #
    # max=10
    # min=-5
    # for i in range(len(genotypes)):
    #     ax[i].get_legend().remove()
    #     ax[i].set_ylim(top=max, bottom=min)
    #     ax[i].set_yticks(np.arange(min, max+1, 5))
    #     # ax.set_title(filename)
    #     ax[i].axhline(y=0, lw=1, ls='--', c='grey')
    #     ax[i].set_ylabel('Change in saccades')
    #     ax[i].set_title(time)
    # fig.savefig(os.path.join(Dir_output, 'saccade_number_change_per_fly_'+filename+time+'.png'))
    # fig.savefig(os.path.join(Dir_output, 'saccade_number_change_per_fly_'+filename+time+'.pdf'), transparent=True, dpi=600)
    # fig.savefig(os.path.join(Dir_output, 'saccade_number_change_per_fly_' + filename + time + '.svg'), transparent=True, dpi=600)
    return 1


def plot_error(genotypes, data_type='', norm=False):
    ## saccadic, smooth response and error comparison
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    fig, ax = plt.subplots(figsize=(3, 3))
    df = pd.DataFrame()

    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    palette = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0) for filename in genotypes]
    # palette = [(0, 0, 0), (1, 0.5, 0.5)]

    for k, Dir in enumerate(genotype_Dir):
        try:
            data = np.load(os.path.join(Dir, genotypes[k]+'_hemi_10_100'+'_cumulative_error_list.npy'), allow_pickle=True)[()]['cumulative_error_list']
            smooth_data = np.load(os.path.join(Dir, genotypes[k]+'_hemi_10_100' + '_smooth_cumulative_error_list.npy'), allow_pickle=True)[()]['cumulative_error_list']
            saccadic_data = np.load(os.path.join(Dir, genotypes[k]+'_hemi_10_100' + '_saccadic_cumulative_error_list.npy'), allow_pickle=True)[()]['cumulative_error_list']
        except:
            data = np.load(os.path.join(Dir, genotypes[k]+'_18_hemi_10_100'+'_cumulative_error_list.npy'), allow_pickle=True)[()]['cumulative_error_list']
            smooth_data = np.load(os.path.join(Dir, genotypes[k]+'_18_hemi_10_100' + '_smooth_cumulative_error_list.npy'), allow_pickle=True)[()]['cumulative_error_list']
            saccadic_data = np.load(os.path.join(Dir, genotypes[k]+'_18_hemi_10_100' + '_saccadic_cumulative_error_list.npy'), allow_pickle=True)[()]['cumulative_error_list']
        for i in range(len(data)):
            if norm:
                new_df = pd.DataFrame.from_dict({'genotype': [genotypes[k]], 'early_error': data[i][0], 'late_error': data[i][1], 'total_error': data[i][2] / data[i][11],
                        'early_FtB': data[i][3], 'late_FtB': data[i][4], 'total_FtB': data[i][5], 'early_BtF': data[i][6], 'late_BtF': data[i][7],'total_BtF': data[i][8],
                        'early_fullrotation': data[i][9], 'late_fullrotation': data[i][10], 'total_fullrotation': data[i][11],'saccadic_fullrotation': saccadic_data[i][11],
                        'smooth_fullrotation': smooth_data[i][11], 'smooth_error': smooth_data[i][2] / smooth_data[i][11],'saccadic_error': saccadic_data[i][2] / saccadic_data[i][11],
                        'smooth_FtB': smooth_data[i][5], 'saccadic_FtB': saccadic_data[i][5], 'smooth_BtF': smooth_data[i][8],'saccadic_BtF': saccadic_data[i][8]})
            else:
                new_df = pd.DataFrame.from_dict({'genotype': [genotypes[k]], 'early_error': data[i][0], 'late_error': data[i][1], 'total_error':data[i][2],
                        'early_FtB': data[i][3], 'late_FtB': data[i][4], 'total_FtB': data[i][5], 'early_BtF': data[i][6], 'late_BtF': data[i][7], 'total_BtF': data[i][8],
                        'early_fullrotation': data[i][9], 'late_fullrotation': data[i][10], 'total_fullrotation': data[i][11], 'saccadic_fullrotation':saccadic_data[i][11],
                        'smooth_fullrotation':smooth_data[i][11], 'smooth_error': smooth_data[i][2], 'saccadic_error': saccadic_data[i][2], 'smooth_FtB':smooth_data[i][5],
                        'saccadic_FtB':saccadic_data[i][5], 'smooth_BtF': smooth_data[i][8], 'saccadic_BtF': saccadic_data[i][8]})
            df = pd.concat([df, new_df], ignore_index=True)

    df = pd.melt(df, id_vars='genotype', var_name='period', value_name='difference')
    # if data_type == 'fullrotation':
    #     sns.stripplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_'+data_type], data=df, ax = ax, palette=palette, size=10, jitter=0.05, alpha=0.5, edgecolor='k', linewidth=0, dodge=True)
    #     annotator = Annotator(ax, data=df, y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_'+data_type], pairs=[(('total_'+data_type, genotypes[0]), ('total_'+data_type, genotypes[1]))])
    annotator = None
    names = ''
    for gen in genotypes:
        names+=gen
    if data_type == 'all':
        sns.pointplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_fullrotation', 'total_BtF', 'total_FtB'], data=df, dodge=0.4, ax=ax,
                      palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_fullrotation', 'total_BtF', 'total_FtB'], data=df, ax = ax, palette=palette, size=7,
                      jitter=0.15, linewidth=0.0, dodge=True)
        pairs = []
        for x in ['total_fullrotation', 'total_BtF', 'total_FtB']:
            for i in range(len(genotypes)-1):
                pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
        annotator = Annotator(ax, data=df, y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_fullrotation', 'total_BtF', 'total_FtB'], pairs=pairs)

    if data_type == 'allsmooth':
        sns.pointplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['smooth_fullrotation', 'smooth_BtF', 'smooth_FtB'], data=df, dodge=0.4, ax=ax,
                      palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['smooth_fullrotation', 'smooth_BtF', 'smooth_FtB'], data=df, ax = ax, palette=palette, size=7,
                      jitter=0.2, linewidth=0.0, dodge=True)
        pairs = []
        for x in ['smooth_fullrotation', 'smooth_BtF', 'smooth_FtB']:
            for i in range(len(genotypes)-1):
                pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
        annotator = Annotator(ax, data=df, y="difference", x='period', hue='genotype', hue_order=genotypes, order=['smooth_fullrotation', 'smooth_BtF', 'smooth_FtB'], pairs=pairs)

    if data_type == 'allsaccadic':
        sns.pointplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_fullrotation', 'total_BtF', 'total_FtB'], data=df, dodge=0.4, ax=ax,
                      palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['saccadic_fullrotation', 'saccadic_BtF', 'saccadic_FtB'], data=df, ax = ax, palette=palette,
                      size=7, jitter=0.15, linewidth=0.0, dodge=True)
        pairs = []
        for x in ['saccadic_fullrotation', 'saccadic_BtF', 'saccadic_FtB']:
            for i in range(len(genotypes)-1):
                pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
        annotator = Annotator(ax, data=df, y="difference", x='period', hue='genotype', hue_order=genotypes, order=['saccadic_fullrotation', 'saccadic_BtF', 'saccadic_FtB'], pairs=pairs)

    if data_type == 'FtBandBtF':
        sns.pointplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_BtF', 'total_FtB'], data=df, dodge=0.4, ax=ax,
                      palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_BtF', 'total_FtB'], data=df, ax = ax, palette=palette,
                      size=7, jitter=0.05, linewidth=0, dodge=True)
        pairs = []
        for x in ['total_BtF', 'total_FtB']:
            for i in range(len(genotypes)-1):
                pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
        annotator = Annotator(ax, data=df, y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_BtF', 'total_FtB'], pairs=pairs)

    if data_type == 'error':
        sns.pointplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_error'], data=df, dodge=0.4, ax=ax,
                      palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(y="difference", x='period', hue='genotype', data=df, ax = ax, palette=palette, size=7, jitter=0.05, linewidth=0,
                      order=['total_error'], hue_order=genotypes, dodge=0.5)
        pairs = []
        for x in ['total_BtF', 'total_FtB']:
            for i in range(len(genotypes)-1):
                pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
        annotator = Annotator(ax, data=df, y="difference", x='period', hue='genotype', order=['total_error'], hue_order=genotypes, pairs=pairs)

    if data_type == 'all_error':
        if len(genotypes)!=1:
            sns.pointplot(y="difference", x='period', hue='genotype', hue_order=genotypes, order=['total_error', 'saccadic_error', 'smooth_error'], data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="difference", x='period', hue='genotype', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=True,
                          order=['total_error', 'saccadic_error', 'smooth_error'], hue_order=genotypes)
            pairs = []
            for x in ['total_error', 'saccadic_error', 'smooth_error']:
                for i in range(len(genotypes) - 1):
                    pairs.append(((x, genotypes[len(genotypes) - 1]), (x, genotypes[i])))
            print(pairs)
            annotator = Annotator(ax, data=df, y="difference", x='period', order=['total_error', 'saccadic_error', 'smooth_error'], pairs=pairs, hue='genotype',hue_order=genotypes)
        else:
            sns.pointplot(y="difference", x='period', order=['total_error', 'saccadic_error', 'smooth_error'], data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="difference", x='period', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=True,
                          order=['total_error', 'saccadic_error', 'smooth_error'])
            t_stat_total, p_value_total = ttest_1samp(df[df['period']=='total_error']['difference'], popmean=0.0)
            print('total',t_stat_total, p_value_total)
            t_stat_saccadic, p_value_saccadic = ttest_1samp(df[df['period'] == 'saccadic_error']['difference'], popmean=0.0)
            print('saccadic', t_stat_saccadic, p_value_saccadic)
            t_stat_smooth, p_value_smooth = ttest_1samp(df[df['period'] == 'smooth_error']['difference'], popmean=0.0)
            print('smooth', t_stat_smooth, p_value_smooth)

    if annotator:
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, names+data_type+'_comparison_all')
    ax.set_title(data_type)
    # ax.legend_.remove()
    # max=3
    # min=-0.5
    max=150
    min=-75
    ax.set_ylim(top=max, bottom=min)
    ax.spines['left'].set_bounds(high=max, low=min)
    ax.spines['bottom'].set_visible(False)
    # ax.set_yticks(np.arange(min, max+0.1, abs(min)))
    ax.set_yticks([-75, -50, 0, 50, 100, 150])
    ax.axhline(y=0, lw=1, ls='dashed', c='grey', alpha=0.5)
    if norm:
        fig.savefig(os.path.join(Dir_output, names+data_type+'_normed_comparison_all.png'), dpi=600)
        fig.savefig(os.path.join(Dir_output, names+data_type+'_normed_comparison_all.pdf'), dpi=600, transparent=True)
        fig.savefig(os.path.join(Dir_output, names + data_type + '_normed_comparison_all.svg'), dpi=600, transparent=True)
    else:
        fig.savefig(os.path.join(Dir_output, names+data_type+'_comparison_all.png'), dpi=600)
        fig.savefig(os.path.join(Dir_output, names+data_type+'_comparison_all.pdf'), dpi=600, transparent=True)
        fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_all.svg'), dpi=600, transparent=True)
    return 1

## plot prediction error line plot
def plot_prediction_error(genotypes):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    fig, ax = plt.subplots(figsize=(3, 3))

    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    palette = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0) for filename in genotypes]
    # palette = [(0, 0, 0), (1, 0.5, 0.5)]

    for k, Dir in enumerate(genotype_Dir):
        for files in os.listdir(Dir):
            if files.endswith('cumulative_error_list.npy') and 'smooth' not in files and 'saccad' not in files:
                print(files)
                data = np.load(os.path.join(Dir, files), allow_pickle=True)[()]['prediction_error']
                other_data = np.load(os.path.join(Dir, files), allow_pickle=True)[()]['cumulative_error_list']
                mean_full_rotation_response=1
                full_rotation_response = 0
                for i in range(len(other_data)):
                    full_rotation_response+=other_data[i][11]
                mean_full_rotation_response=full_rotation_response/len(other_data)
                ax.plot(savgol_filter(data, 31, 2), c=palette[k], lw=1.5, label=genotypes[k])

    ax.axhline(y=0, lw=1, ls='dashed', c='grey', alpha=0.5)
    ax.axvline(x=120, lw=1, ls='dashed', c='grey', alpha=0.5)
    # ax.axvline(x=150, lw=1, ls='dashed', c='grey', alpha=0.5)
    # ax.axvline(x=180, lw=1, ls='dashed', c='grey', alpha=0.5)
    # ax.axvline(x=120+180, lw=1, ls='dashed', c='grey', alpha=0.5)
    # max = 1.5
    # min = -0.5
    max = 150
    min = -100
    ax.set_ylim(top=max, bottom=min)
    ax.set_xticks([0, 120, 420])
    ax.set_xticklabels(['-2', '0', '5'])
    ax.spines['bottom'].set_bounds(high=420, low=0)
    ax.spines['left'].set_bounds(high=150, low=min)
    ax.set_yticks([-100, -50, 0, 50, 100, 150])
    ax.set_yticks(np.arange(min, max+0.01, 50))
    ax.legend()
    names=''
    for gen in genotypes:
        names += gen
    ax.set_ylabel('Prediction Error [Full - (FtB + BtF)]/Full')
    fig.savefig(os.path.join(Dir_output, names+'prediction_error_norm.png'), dpi=600)
    fig.savefig(os.path.join(Dir_output, names+'prediction_error_norm.pdf'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir_output, names + 'prediction_error_norm.svg'), dpi=600, transparent=True)
    return 1


def plot_summary_scatter(genotypes):
    # # scatter plot summarizing everything
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'

    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    color_list = [[genotype_colors[filename]]*3 if filename in genotype_colors else [(0, 0, 0)]*3 for filename in genotypes]
    # color_list =  [[(0, 0, 0)]*3, [(1, 0.5, 0.5)]*3, [(0.8, 0.5, 1)]*3]

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    marker_list = ['<', 'D', 'o']

    facecolor_list = ['none', 'none', 'none']
    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    for k, Dir in enumerate(genotype_Dir):
        print(Dir)
        for files in os.listdir(Dir):
            if 'saccadic_smooth_response' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]

        saccadic_angular_speed_eachfly_list = data['saccadic_response_eachfly']
        smooth_angular_speed_eachfly_list = data['smooth_response_eachfly']
        saccadic_angular_speed_sem_eachfly_list = data['saccadic_response_sem_eachfly']
        smooth_angular_speed_sem_eachfly_list = data['smooth_response_sem_eachfly']

        for i in range(len(saccadic_angular_speed_eachfly_list)):
            ## eacf fly
            ax.scatter(x=saccadic_angular_speed_eachfly_list[i], y=smooth_angular_speed_eachfly_list[i], marker=marker_list[i], s=20, linewidths=1.5, edgecolors=color_list[k][i],
                       facecolors=facecolor_list[i], label=labels[i], alpha=0.4, zorder=2)
            ## average of flies
            ax.scatter(x=np.mean(saccadic_angular_speed_eachfly_list[i]), y=np.mean(smooth_angular_speed_eachfly_list[i]), marker=marker_list[i], s=70, linewidths=0, edgecolors=color_list[k][-1],
                       facecolors=color_list[k][-1], alpha=0.85, zorder=2.5)
            ## global average and sem
            ax.errorbar(x=np.mean(saccadic_angular_speed_eachfly_list[i]), y=np.mean(smooth_angular_speed_eachfly_list[i]), xerr=sem(saccadic_angular_speed_eachfly_list[i]),
                        yerr=sem(smooth_angular_speed_eachfly_list[i]), elinewidth=1, ecolor=color_list[k][-1], zorder=3, alpha=0.5, capsize=1)
            ## per fly average and sem
            ax.errorbar(x=saccadic_angular_speed_eachfly_list[i], y=smooth_angular_speed_eachfly_list[i], xerr=saccadic_angular_speed_sem_eachfly_list[i],
                        yerr=smooth_angular_speed_sem_eachfly_list[i], elinewidth=0.5, ecolor=color_list[k][-1], zorder=3, alpha=0.3, capsize=0.5, fmt='none')

    max = 150
    # ax.plot([-max, max], [-max, max], ls=':', lw=0.5, c='grey'
    ax.axvline(x=0, ls='--', c='grey', lw=0.5, alpha=0.4)
    ax.axhline(y=0, ls='--', c='grey', lw=0.5, alpha=0.4)
    ax.spines['bottom'].set_position(('outward', 2))
    ax.spines['left'].set_position(('outward', 2))
    ax.set_xlabel('Saccadic Response')
    ax.set_ylabel('Smooth Response')
    ax.set_ylim(top=max, bottom=-75)
    ax.set_xlim(left=-75, right=max)
    ax.set_yticks(np.arange(-75, max+1, 25))
    ax.set_xticks(np.arange(-75, max+1, 25))
    # ax.legend()
    # ax.set_title(genotypes[0] + '_' +genotypes[1])
    ax.set_title(genotypes[0])
    names = ''
    for gen in genotypes:
        names+=gen
    ax.set_aspect('equal')
    fig.savefig(os.path.join(Dir_output, names + '_summary.png'))
    fig.savefig(os.path.join(Dir_output, names + '_summary.pdf'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir_output, names + '_summary.svg'), dpi=600, transparent=True)
    return 1


def plot_saccade_summary_scatter(genotypes):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    # scatter plot summarizing anti and syn saccades
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    color_list = [[genotype_colors[filename]]*3 if filename in genotype_colors else [(0, 0, 0)]*3 for filename in genotypes]
    # color_list = [[(0, 0, 0)]*3, [(1, 0.5, 0.5)]*3, [(0.8, 0.5, 1)]*3]

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    marker_list = ['<', 'D', 'o']

    facecolor_list = ['none', 'none', 'none']
    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    for k, Dir in enumerate(genotype_Dir):
        print(Dir)
        for files in os.listdir(Dir):
            if 'saccadic_smooth_response' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]
        anti_saccadic_turns = data['anti_saccadic_turns_eachfly']
        syn_saccadic_turns = data['syn_saccadic_turns_eachfly']
        anti_saccadic_turns_sem_eachfly = data['anti_saccadic_sem_turns_eachfly']
        syn_saccadic_turns_sem_eachfly = data['syn_saccadic_sem_turns_eachfly']

        for i in range(len(anti_saccadic_turns)):
            print(anti_saccadic_turns_sem_eachfly[i])
            ax.scatter(x=syn_saccadic_turns[i], y=anti_saccadic_turns[i], marker=marker_list[i], s=20, linewidths=1.5, edgecolors=color_list[k][i],
                       facecolors=facecolor_list[i], label=labels[i], alpha=0.4, zorder=2)
            ax.scatter(x=np.mean(syn_saccadic_turns[i]), y=np.mean(anti_saccadic_turns[i]), marker=marker_list[i], s=70, linewidths=0, edgecolors=color_list[k][-1],
                       facecolors=color_list[k][-1], alpha=0.85, zorder=2.5)
            ax.errorbar(x=np.mean(syn_saccadic_turns[i]), y=np.mean(anti_saccadic_turns[i]), xerr=sem(syn_saccadic_turns[i]), yerr=sem(anti_saccadic_turns[i]),
                        elinewidth=1, ecolor=color_list[k][-1], zorder=3, alpha=0.5, capsize=1)
            ax.errorbar(x=syn_saccadic_turns[i], y=anti_saccadic_turns[i], xerr=syn_saccadic_turns_sem_eachfly[i], yerr=anti_saccadic_turns_sem_eachfly[i],
                        elinewidth=0.5, ecolor=color_list[k][-1], zorder=3, alpha=0.3, capsize=0.5, fmt='none')

    max = 200
    ax.plot([0, 400], [0, -400], c='grey', lw=0.5, alpha=0.4, ls='--')
    # ax.axvline(x=0, ls='--', c='grey', lw=0.5, alpha=0.4)
    # ax.axhline(y=0, ls='--', c='grey', lw=0.5, alpha=0.4)
    ax.xaxis.tick_top()
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(True)
    ax.spines['top'].set_position(('outward', 2))
    ax.spines['left'].set_position(('outward', 2))
    ax.set_xlabel('Syn-saccadic response')
    ax.set_ylabel('Anti-saccadic response')
    ax.set_ylim(top=-25, bottom=-400)
    ax.set_xlim(left=-10, right=400)
    ax.set_yticks(np.arange(-400, 1, 100))
    ax.set_xticks(np.arange(0, 401, 100))
    ax.spines['top'].set_bounds(low=0, high=400)
    ax.spines['left'].set_bounds(low=-400, high=0)
    # ax.legend()
    # ax.set_title(genotypes[0] + '_' +genotypes[1])
    names = ''
    for gen in genotypes:
        names+=gen
    ax.set_aspect('equal')
    fig.savefig(os.path.join(Dir_output, names + 'saccade_summary.png'))
    fig.savefig(os.path.join(Dir_output, names + 'saccade_summary.pdf'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir_output, names + 'saccade_summary.svg'), dpi=600, transparent=True)
    return 1

def plot_saccade_shapes(genotypes):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    # scatter plot summarizing anti and syn saccades
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    color_list = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0, 0) for filename in genotypes]

    fig, ax = plt.subplots(2, 2, figsize=(5, 5), sharex=True)

    facecolor_list = ['none', 'none', 'none']
    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    for k, Dir in enumerate(genotype_Dir):
        for files in os.listdir(Dir):
            if 'saccade_shape_hemi' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]

        ax[1, 1].plot(np.mean(np.vstack((-1*np.array(data['ccw_syn_angular_speed']), 1*np.array(data['cw_syn_angular_speed']))), axis=0), lw=3, color=color_list[k], alpha=0.6)
        ax[0, 1].plot(np.mean(np.vstack((-1*np.array(data['ccw_syn_opto_response']), 1*np.array(data['cw_syn_opto_response']))), axis=0), lw=3, color=color_list[k], alpha=0.6)
        ax[1, 0].plot(np.mean(np.vstack((1*np.array(data['ccw_anti_angular_speed']), -1*np.array(data['cw_anti_angular_speed']))), axis=0), lw=3, color=color_list[k], alpha=0.6)
        ax[0, 0].plot(np.mean(np.vstack((1*np.array(data['ccw_anti_opto_response']), -1*np.array(data['cw_anti_opto_response']))), axis=0), lw=3, color=color_list[k], alpha=0.6)

        ax[1, 0].set_ylim(top=50, bottom=-50)
        ax[1, 1].set_ylim(top=50, bottom=-50)

        ax[0, 1].set_ylim(top=600, bottom=-50)
        ax[0, 0].set_ylim(top=600, bottom=-50)

        ax[1, 0].set_xticks([0, 5, 10, 15, 20])
        ax[1, 0].set_xticklabels(['-160', '-80', '0', '80', '160'])

    names = ''
    for gen in genotypes:
        names+=gen
    fig.savefig(os.path.join(Dir_output, names + 'saccade_shape.png'))
    fig.savefig(os.path.join(Dir_output, names + 'saccade_shape.pdf'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir_output, names + 'saccade_shape.svg'), dpi=600, transparent=True)
    return 1


def plot_straightness_per_fly(genotypes):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures\NewPlotsSuppl2024'

    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    genotype_Dir = [r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\CantonS_18\Hemisphere\Contrast_100\ProcessedData']
    color_list = [np.append(genotype_colors[filename], 0.5) for filename in genotypes]
    print(color_list)
    fig, ax = plt.subplots(figsize=(2.5*(9/6), 2.5*(9/6)))
    labels = ['pre', 'FullRotation', 'BackToFront', 'FrontToBack']
    df = []
    for k, Dir in enumerate(genotype_Dir):
        print(Dir)
        for files in os.listdir(Dir):
            if '_straightness.csv' in files:
                print(files)
                path = os.path.join(Dir, files)
                break
        data = pd.read_csv(path)
        df.append(data)
    df = pd.concat(df, ignore_index=True)
    df.to_csv(os.path.join(Dir_output, 'output.csv'))
    names = ''
    for gen in genotypes:
        names+=gen

    if len(genotypes) == 1:
        sns.pointplot(data=df, y='straightness', x='stim_type', order=labels, dodge=0.4, ax=ax, palette=color_list, estimator=np.mean,
                      errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(data=df, y='straightness', x='stim_type', order=labels, palette=color_list, size=7, jitter=0.15, linewidth=0.0,
                      dodge=False, ax=ax)
        pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack', 'pre'], 2))
        annotator = Annotator(ax, data=df, y='straightness', x='stim_type', order=labels, pairs=pairs)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, names + '_straightness')
    else:
        sns.pointplot(data=df, y='straightness', x='stim_type', order=labels, hue='genotype', hue_order=genotypes, dodge=0.4, ax=ax, palette=color_list, estimator=np.mean,
                      errorbar='se', join=False, markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
        sns.stripplot(data=df, y='straightness', x='stim_type', order=labels, hue='genotype', hue_order=genotypes, palette=color_list, size=7, jitter=0.15, linewidth=0.0,
                      dodge=True, ax=ax)
        pairs = []
        for x in ['FrontToBack', 'BackToFront', 'FullRotation', 'pre']:
            for i in range(1, len(genotypes)):
                pairs.append(((x, genotypes[0]), (x, genotypes[i])))
        annotator = Annotator(ax, data=df, y='straightness', x='stim_type', hue='genotype', hue_order=genotypes, order=labels, pairs=pairs)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, names + '_straightness')

    max=1
    min=0.5
    ax.set_ylim(top=max, bottom=min)
    ax.spines['left'].set_bounds(high=max, low=min)
    ax.spines['bottom'].set_visible(False)
    ax.set_yticks(np.arange(min, max+0.01, 0.1))

    fig.savefig(os.path.join(Dir_output, names + '_straightness.png'))
    fig.savefig(os.path.join(Dir_output, names + '_straightness.pdf'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir_output, names + '_straightness.svg'), dpi=600, transparent=True)
    return 1


def plot_straightness(genotypes):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    df = pd.DataFrame({'straightness':[0], 'genotype':['a'], 'stim_type':['a']})

    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    color_list = [np.append(genotype_colors[filename], 0.5) for filename in genotypes]
    # color_list = [[0.5, 0.5, 0.5, 0.5] for filename in genotypes]

    fig, ax = plt.subplots(figsize=(2.5, 2.5))

    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    for k, Dir in enumerate(genotype_Dir):
        for files in os.listdir(Dir):
            if '_straightness' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]
        print(len(data['prestim']))

        default_straightness = []
        for i in range(len(data['stim'])):
            df1 = pd.DataFrame({'straightness': data['stim'][i], 'genotype': [genotypes[k]] * data['stim'][i].shape[0], 'stim_type': [labels[i]] * data['stim'][i].shape[0]})
            df = pd.concat([df, df1], ignore_index=True)
            default_straightness+=list(data['prestim'][i])

        df1 = pd.DataFrame({'straightness': default_straightness, 'genotype': [genotypes[k]] * len(default_straightness), 'stim_type': ['prestim'] * len(default_straightness)})
        df = pd.concat([df, df1], ignore_index=True)
    df = df.drop([0])

    labels = ['FullRotation', 'BackToFront', 'FrontToBack']
    if len(genotypes)==1:
        sns.violinplot(data=df, y='straightness', x='stim_type', order=labels, hue='genotype', hue_order=genotypes, palette=color_list, linewidth=0.5, scale='width', cut=0, inner=None)
        sns.boxplot(data=df, y='straightness', x='stim_type', order=labels, linewidth=1, color='grey', showfliers=False, width=0.2, showcaps=False,
                    boxprops={"facecolor": (.0, .0, .0, .5), 'zorder': 2, 'linewidth':0.0},
                    medianprops={"color": "white"})
    else:
        sns.violinplot(data=df, y='straightness', x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack', 'prestim'], hue='genotype', hue_order=genotypes, palette=color_list, linewidth=0.5, scale='width', cut=0,
                       inner=None, split=True)
    # ax.get_legend().remove()
    # ax.set_ylim(top=1, bottom=0.25)
    names = ''
    for gen in genotypes:
        names += gen

    fig.savefig(os.path.join(Dir_output, names + '_straightness.png'), dpi=600)
    fig.savefig(os.path.join(Dir_output, names + '_straightness.pdf'), dpi=600, transparent=True)
    fig.savefig(os.path.join(Dir_output, names + '_straightness.svg'), dpi=600, transparent=True)
    return 1


def plot_syn_anti_saccade_alltrials(genotypes, data_type=''):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    palette = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0) for filename in genotypes]
    # palette = [(0, 0, 0), (1, 0.5, 0.5)]

    fig, ax = plt.subplots(figsize=(2.5, 2.5))

    df = pd.DataFrame({'anti_saccadic_turns': [0], 'syn_saccadic_turns': [0], 'smooth_turns': [0], 'saccade_turns': [0],'genotype': ['control'], 'stim_type':['red']})
    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    print(genotype_Dir)
    for k, Dir in enumerate(genotype_Dir):
        for files in os.listdir(Dir):
            if 'saccadic_smooth_response' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]

        anti_saccadic_turns = data['anti_saccadic_turns_all']
        syn_saccadic_turns = data['syn_saccadic_turns_all']
        smooth_turns = data['smooth_turns_all']
        saccade_turns = data['saccadic_turns_all']

        for i in range(3):
            df1 = pd.DataFrame({'anti_saccadic_turns': anti_saccadic_turns[i], 'syn_saccadic_turns': syn_saccadic_turns[i], 'smooth_turns': smooth_turns[i], 'saccade_turns': saccade_turns[i],
                                'genotype':[genotypes[k]]*len(anti_saccadic_turns[i]), 'stim_type': [labels[i]]*len(anti_saccadic_turns[i])})
            df = pd.concat([df, df1], ignore_index=True)
    annotator = None
    # print(df)
    # print(df.loc[df['stim_type']=='FrontToBack']['anti_saccadic_turns'])
    # print(df.loc[df['stim_type']=='BackToFront']['syn_saccadic_turns'])
    # print(mannwhitneyu(np.abs(df.loc[df['stim_type'] == 'FrontToBack']['anti_saccadic_turns']), df.loc[df['stim_type'] == 'BackToFront']['syn_saccadic_turns']))
    # return 1
    if data_type == 'all_anti_saccadic':
        if len(genotypes)==1:
            sns.pointplot(y="anti_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y="anti_saccadic_turns", x='stim_type', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="anti_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y="anti_saccadic_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            annotator = Annotator(ax, data=df, y="anti_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    if data_type == 'all_syn_saccadic':
        if len(genotypes)==1:
            sns.pointplot(y="syn_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y="syn_saccadic_turns", x='stim_type', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="syn_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="syn_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y="syn_saccadic_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            annotator = Annotator(ax, data=df, y="syn_saccadic_turns", x='stim_type',order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    if data_type == 'all_smooth':
        if len(genotypes)==1:
            sns.pointplot(y="smooth_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y='smooth_turns', x='stim_type', data=df, ax=ax, palette=palette, size=0.75,  linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="smooth_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="smooth_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y='smooth_turns', x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            annotator = Annotator(ax, data=df, y="smooth_turns", x='stim_type',order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    if data_type == 'all_saccadic':
        if len(genotypes)==1:
            sns.pointplot(y="saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y="saccade_turns", x='stim_type', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="anti_saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.swarmplot(y="saccade_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=0.75, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            annotator = Annotator(ax, data=df, y="saccade_turns", x='stim_type',order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    names = ''
    for gen in genotypes:
        names += gen

    ax.set_title(data_type)
    # ax.legend_.remove()
    # max=3
    # min=-0.5
    if data_type == 'all_syn_saccadic':
        max = 500
        min = -25
        ax.set_yticks([0, 100, 200, 300, 400, 500])
    elif data_type == 'all_anti_saccadic':
        max = 25
        min = -500
        ax.set_yticks([-500, -400, -300, -200, -100, 0])
    elif data_type == 'all_smooth':
        max = 600
        min = -300
        ax.set_yticks([-300, -150, 0, 250, 500, 600])
    ax.set_ylim(top=max, bottom=min)
    ax.spines['left'].set_bounds(high=max, low=min)
    ax.spines['bottom'].set_visible(False)
    ax.legend([], [])
    # ax.set_yticks(np.arange(min, max+0.1, abs(min)))

    if annotator:
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, names + data_type + '_comparison_alltrials')

    ax.axhline(y=0, lw=1, ls='dashed', c='grey', alpha=0.5)
    fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_alltrials.png'), dpi=600)
    fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_alltrials.pdf'), dpi=600, transparent=True)
    # fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_alltrials.svg'), dpi=600, transparent=True)
    return 1


def plot_syn_anti_saccade(genotypes, data_type=''):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures\NewPlotsSuppl2024'
    print(genotype_folders)
    print(genotypes)
    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    print(genotype_Dir)
    palette = [genotype_colors[filename] if filename in genotype_colors else (0, 0, 0) for filename in genotypes]
    cmap = make_color_pastel(ListedColormap(all_color_data['PiYG5']))
    # palette = [cmap(0.95), cmap(0.05)]
    # palette = [(0, 0, 0), (1, 0.5, 0.5)]

    fig, ax = plt.subplots(figsize=(2.5, 2.5))

    df = pd.DataFrame({'anti_saccadic_turns': [0], 'syn_saccadic_turns': [0], 'smooth_turns': [0], 'saccade_turns': [0],'genotype': ['control'], 'stim_type':['red']})
    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    print(genotype_Dir)
    for k, Dir in enumerate(genotype_Dir):
        for files in os.listdir(Dir):
            if 'saccadic_smooth_response' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]
        anti_saccadic_turns = data['anti_saccadic_turns_eachfly']
        syn_saccadic_turns = data['syn_saccadic_turns_eachfly']
        # smooth_turns = data['smooth_response_eachfly']
        # saccade_turns = data['saccadic_response_eachfly']
        # changed on 27/11
        smooth_turns = data['smooth_turns_eachfly']
        saccade_turns = data['saccadic_turns_eachfly']

        anti_saccadic_turns_all = data['anti_saccadic_turns_eachfly']
        syn_saccadic_turns_all = data['syn_saccadic_turns_eachfly']
        smooth_turns_all = data['smooth_turns_eachfly']
        saccade_turns_all = data['saccadic_turns_eachfly']

        for i in range(3):
            df1 = pd.DataFrame({'anti_saccadic_turns': anti_saccadic_turns[i], 'syn_saccadic_turns': syn_saccadic_turns[i], 'smooth_turns': smooth_turns[i], 'saccade_turns': saccade_turns[i],
                                'genotype':[genotypes[k]]*len(anti_saccadic_turns[i]), 'stim_type': [labels[i]]*len(anti_saccadic_turns[i])})
            df = pd.concat([df, df1], ignore_index=True)
    annotator = None
    # print(df)
    # print(df.loc[df['stim_type']=='FrontToBack']['anti_saccadic_turns'])
    # print(df.loc[df['stim_type']=='BackToFront']['syn_saccadic_turns'])
    # print(mannwhitneyu(np.abs(df.loc[df['stim_type'] == 'FrontToBack']['anti_saccadic_turns']), df.loc[df['stim_type'] == 'BackToFront']['syn_saccadic_turns']))
    # return 1
    if data_type == 'all_anti_saccadic':
        if len(genotypes)==1:
            sns.pointplot(y="anti_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="anti_saccadic_turns", x='stim_type', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="anti_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=[(0.2, 0.2, 0.2) for i in range(len(genotypes))], estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="anti_saccadic_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            for y in genotypes:
                for label in list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2)):
                    pairs.append(((label[0], y), (label[1], y)))
            annotator = Annotator(ax, data=df, y="anti_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    if data_type == 'all_syn_saccadic':
        if len(genotypes)==1:
            sns.pointplot(y="syn_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="syn_saccadic_turns", x='stim_type', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="syn_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="syn_saccadic_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=[(0.2, 0.2, 0.2) for i in range(len(genotypes))], estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="syn_saccadic_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            for y in genotypes:
                for label in list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2)):
                    pairs.append(((label[0], y), (label[1], y)))
            annotator = Annotator(ax, data=df, y="syn_saccadic_turns", x='stim_type',order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    if data_type == 'all_smooth':
        if len(genotypes)==1:
            sns.pointplot(y="smooth_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y='smooth_turns', x='stim_type', data=df, ax=ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="smooth_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="smooth_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=[(0.2, 0.2, 0.2) for i in range(len(genotypes))], estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y='smooth_turns', x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            for y in genotypes:
                for label in list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2)):
                    pairs.append(((label[0], y), (label[1], y)))
            annotator = Annotator(ax, data=df, y="smooth_turns", x='stim_type',order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    if data_type == 'all_saccadic':
        if len(genotypes)==1:
            sns.pointplot(y="saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], data=df, dodge=0.4, ax=ax,
                          palette=palette, estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="saccade_turns", x='stim_type', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=False,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'])
            pairs = list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2))
            annotator = Annotator(ax, data=df, y="anti_saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], pairs=pairs)
        else:
            sns.pointplot(y="saccade_turns", x='stim_type', order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', data=df, dodge=0.4, ax=ax,
                          palette=[(0.2, 0.2, 0.2) for i in range(len(genotypes))], estimator=np.mean, errorbar='se', join=False, markers="_", errwidth=1, scale=2, hue_order=genotypes)
            plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
            plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
            sns.stripplot(y="saccade_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=palette, size=7, jitter=0.15, linewidth=0.0, dodge=True,
                          order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)
            pairs = []
            for x in ['FullRotation', 'BackToFront', 'FrontToBack']:
                for i in range(len(genotypes)-1):
                    pairs.append(((x, genotypes[len(genotypes)-1]), (x, genotypes[i])))
            for y in genotypes:
                for label in list(itertools.combinations(['FullRotation', 'BackToFront', 'FrontToBack'], 2)):
                    pairs.append(((label[0], y), (label[1], y)))
            annotator = Annotator(ax, data=df, y="saccade_turns", x='stim_type',order=['FullRotation', 'BackToFront', 'FrontToBack'], hue='genotype', hue_order=genotypes, pairs=pairs)

    names = ''
    for gen in genotypes:
        names += gen

    ax.set_title(data_type)
    # ax.legend_.remove()
    # max=3
    # min=-0.5
    if data_type == 'all_syn_saccadic':
        max = 400
        min = -25
        ax.set_yticks([0, 100, 200, 300, 400, 500])
    elif data_type == 'all_anti_saccadic':
        max = 25
        min = -400
        ax.set_yticks([-500, -400, -300, -200, -100, 0])
    elif data_type == 'all_smooth':
        max = 200
        min = -100
        ax.set_yticks([-100, 0, 100, 200])
    ax.set_ylim(top=max, bottom=min)
    ax.spines['left'].set_bounds(high=max, low=min)
    ax.spines['bottom'].set_visible(False)
    # ax.set_yticks(np.arange(min, max+0.1, abs(min)))

    if annotator:
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
        _, annotations = annotator.apply_and_annotate()
        save_stats_data(annotations, Dir_output, names + data_type + '_comparison_all')

    ax.axhline(y=0, lw=1, ls='dashed', c='grey', alpha=0.5)
    fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_all.png'), dpi=600)
    fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_all.pdf'), dpi=600, transparent=True)
    # fig.savefig(os.path.join(Dir_output, names + data_type + '_comparison_all.svg'), dpi=600, transparent=True)
    return 1

def plot_all_respnses(genotypes):
    Dir_output = r'D:\Roshan\Experiments\Optomotor\NewDataSmallArena\FinalData\Figures'
    df = pd.DataFrame({'straightness':[0], 'genotype':['a'], 'stim_type':['a']})

    genotype_Dir = [genotype_folders[filename] for filename in genotypes]
    color_list = [np.append(genotype_colors[filename], 0.5) for filename in genotypes]

    df = pd.DataFrame({'opto_turns_response': [0], 'angular_speed_response': [0], 'genotype': ['control'], 'stim_type':['red']})
    labels = ['FrontToBack', 'BackToFront', 'FullRotation']
    for k, Dir in enumerate(genotype_Dir):
        for files in os.listdir(Dir):
            if '_total_data_hemi' in files:
                path = os.path.join(Dir, files)
                break
        data = np.load(path, allow_pickle=True)[()]

        for i in range(3):
            df1 = pd.DataFrame({'opto_turns_response': data['opto_response_list'][i], 'angular_speed_response': data['angular_speed_list'][i],
                                'genotype': [genotypes[k]]*len(data['opto_response_list'][i]), 'stim_type': [labels[k]]*len(data['opto_response_list'][i])})
            df = pd.concat([df, df1], ignore_index=True)
    df = df.dop([0])

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    sns.kdeplot(y="opto_response_turns", x='stim_type', hue='genotype', data=df, ax = ax, palette=color_list, size=7*(1.11/0.92), jitter=0.15, linewidth=0.0, dodge=True,
                      order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    sns.kdeplot(y="angular_speed_response", x='stim_type', hue='genotype', data=df, ax = ax, palette=color_list, size=7*(1.11/0.92), jitter=0.15, linewidth=0.0, dodge=True,
                      order=['FullRotation', 'BackToFront', 'FrontToBack'], hue_order=genotypes)

    return 1