import mat73
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import matplotlib.gridspec as gridspec
import numpy as np
import json
import matplotlib.colors as cm
import pandas as pd
from scipy.stats import sem
from scipy.interpolate import interp1d
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
from scipy.stats import mannwhitneyu, ttest_ind
from matplotlib.transforms import Bbox
from scipy.integrate import trapz
import seaborn as sns
from statannotations.Annotator import Annotator
from scipy.ndimage import zoom

global color_dict
color_dict = json.load(r'path\to\color_dict.json') ## a dictionary of colors for different genotypes
shades_ND = color_dict['FlpND']
shades_D = color_dict['FlpD']
shades_ShakB2 = color_dict['ShakB2']
shades_DB331 = color_dict['FlpNDxDB331']

global Dir, Dir_output
Dir = r'path\to\data'
Dir_output = r'path\to\save_destination'

## helper functions
def sph2cart(az,elev,r):
    z = r * np.sin(elev)
    rcoselev = r * np.cos(elev)
    x = rcoselev * np.cos(az)
    y = rcoselev * np.sin(az)
    return [x, y, z]

def expand_grid(ncol, nrow):
    hor_degree_span_half = 70
    ver_degree_span_half = 342 * hor_degree_span_half / 608
    azvals = np.linspace(-hor_degree_span_half, hor_degree_span_half, ncol)
    elvals = np.linspace(-ver_degree_span_half, ver_degree_span_half, nrow)
    [azmesh, elmesh] = np.meshgrid(azvals, elvals)

    # wrap the grid on the sphere
    # get [x,y] coordinatas in 3D corresponding to the grid of scanlines
    azi = np.pi * azmesh / 180
    eli = np.pi * elmesh / 180
    [z3d, y3d, x3d] = sph2cart(eli, azi, 1)
    spherical_coords_grid_cart = np.stack((x3d, y3d, z3d), axis=2)

    return spherical_coords_grid_cart

def expand_grid_spherical(ncol, nrow):
    hor_degree_span_half = 70
    ver_degree_span_half = 342 * hor_degree_span_half / 608
    azvals = np.linspace(-hor_degree_span_half, hor_degree_span_half, ncol)
    elvals = np.linspace(-ver_degree_span_half, ver_degree_span_half, nrow)
    [azmesh, elmesh] = np.meshgrid(azvals, elvals)

    # # wrap the grid on the sphere
    # # get [x,y] coordinatas in 3D corresponding to the grid of scanlines
    azi = np.pi * azmesh / 180
    eli = np.pi * elmesh / 180
    # [z3d, y3d, x3d] = sph2cart(eli, azi, 1)
    spherical_coords_grid_cart = np.stack((azi, eli), axis=2)

    return spherical_coords_grid_cart

def significance(data):
    stars=[]
    for d in data:
        if d>0.05:
            stars.append(0)
        elif 0.01<d<0.05:
            stars.append(1)
        elif 0.001<d<0.01:
            stars.append(2)
        elif d<0.001:
            stars.append(3)
    return stars

def create_color_gradient(color1, color2, n):
    color_list = []
    t = np.linspace(1, 0, n)
    for i in t:
        color_list.append(np.array(color1) + np.multiply(np.subtract(color2, color1), i))
    return color_list


## plotting functions
## make video of RFs changing from FlpND to FlpD
def make_RF_video():
    ## make RF video
    data_FlpND = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + 'HSE' + '_' + 'FlpND' + '.mat'))
    data_FlpD = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + 'HSE' + '_' + 'FlpD' + '.mat'))
    RF_FlpND = data_FlpND['grp_data']['meanRF']
    RF_FlpD = data_FlpD['grp_data']['meanRF']

    def my_func(x,y, n=60):
        return np.interp(np.linspace(0, 1, n), [0, 1], [x, y])
    new_interp = np.vectorize(my_func,excluded=['n'], signature='(),()->(k)')

    all_transition_RFs = new_interp(RF_FlpND, RF_FlpD)

    grid = data_FlpND['grp_data']['spherical_grid']
    fig, ax = plt.subplots(figsize=(grid.shape[1]/ 5, grid.shape[0]/ 5))
    data = all_transition_RFs
    artists = []
    grid = data_FlpD['grp_data']['deformed_grid']
    x = grid[:, :, 0]
    y = grid[:, :, 1]
    n = 10
    ax.set_xlim(x[0, -1] + (x[0, -1] / n), x[0, 0] + (x[0, 0] / n))
    # y_min = np.amin(y)
    # y_max = np.amax(y)
    # m=9
    # ax_main.set_ylim(y_min+(y_min/m), y_max+(y_max/m))
    # ax_yDist.set_ylim(y_min+(y_min/m), y_max+(y_max/m))
    ##
    ax.set_xticks([x[0, -1], 0, x[0, 0]])
    ax.set_yticks([y[-1, 6], 0, y[0, 6]])

    color_list1 = create_color_gradient([1, 1, 1], color_dict['FlpND'][0], 45)
    color_list2 = create_color_gradient(color_dict['FlpD'][0], [1, 1, 1], 45)
    color_list = color_list1[:30] + color_list2[15:]
    for i in range(90):
        if i < 25:
            k = 0
        elif i>=85:
            k = -1
        else:
            print(i)
            k = i-25
        u = all_transition_RFs[:, :, 0, k]*-1
        v = all_transition_RFs[:, :, 1, k]
        # artists.append((ax.quiver(x, y, u, v, width=0.005, scale=1.5, headwidth=4, headlength=6, headaxislength=5, color=cmap((60-i)/60),
        #                alpha=1, minlength=1, zorder=2.5),))
        artists.append((ax.quiver(x, y, u, v, width=0.005, scale=1.4, headwidth=3, headlength=6, headaxislength=5, color=color_list[k],
                    alpha=1, minlength=1, zorder=2.5),))
        # ax.set_facecolor('k')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

    im_ani = animation.ArtistAnimation(fig, artists, repeat_delay=0, blit=True)
    im_ani.save(os.path.join(Dir_output, 'RF_change.mp4'), writer='ffmpeg', fps=30, dpi=300)
    return 1


## Fig4a,b,c,d
def make_RF_plots(cellname_list=['HSN'], direction_list=['neg_x'], flies_to_plot_list=[['FlpND', 'FlpD']], response=0.15, min_response=-0.15, 
                  draw_contour=False, expansion_factor=5, draw_contour_only_for_horizontal=False, draw_grid_colors=False, contour_levels=[0.6, 1], spiking=False):
    for cellname in cellname_list:
        for direction in direction_list:
            for flies_to_plot in flies_to_plot_list:
                response = 0.15
                min_response = -0.15
                # response = 0.15
                # min_response = 0
                data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname + '_' + 'FlpND' + '.mat'))
                grid = data['grp_data']['spherical_grid']
                fig = plt.figure(figsize=((grid.shape[1] + 2) / 5, (grid.shape[0] + 2) / 5))
                print((grid.shape[1] + 2) / 6, (grid.shape[0] + 2) / 6)
                gs = gridspec.GridSpec(grid.shape[0] + 2, grid.shape[1] + 2, wspace=0, hspace=0)
                ax_main = plt.subplot(gs[2:grid.shape[0] + 2, :grid.shape[1]])
                ax_xDist = plt.subplot(gs[:2, :grid.shape[1]], sharex=ax_main)
                ax_yDist = plt.subplot(gs[2:grid.shape[0] + 2, grid.shape[1]:], sharey=ax_main)

                filename = ''
                for names in flies_to_plot:
                    filename += names
                    if spiking:
                        print('spiking')
                        data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname+'_' + names + "_spike" + '.mat'))
                    else:
                        data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname+'_' + names + '.mat'))

                    grid = data['grp_data']['deformed_grid']
                    if direction == '':
                        RF = data['grp_data']['meanRF']
                        hor = data['grp_data']['RF_hor_profile']
                        vert = data['grp_data']['RF_ver_profile']
                        RF_amp = data['grp_data']['RF_ampmean']
                        cnt_levels = data['grp_data']['cntr_levels']
                    else:
                        RF = data['grp_data']['meanRF'+'_'+direction]
                        hor = data['grp_data']['RF_hor_profile'+'_'+direction]
                        vert = data['grp_data']['RF_ver_profile'+'_'+direction]
                        hor_real = data['grp_data']['RF_hor_vec_profile'+'_'+direction][:, 0]*-1 #only take RFx means
                        vert_real = data['grp_data']['RF_ver_vec_profile'+'_'+direction][:, 0]*-1
                        hor = hor_real
                        vert = vert_real
                        RF_amp = data['grp_data']['RF_ampmean'+'_'+direction]
                        # cnt_levels = data['grp_data']['cntr_levels'+'_'+direction]
                    print(RF.shape)
                    h_span = data['grp_data']['hspan_half']
                    v_span = data['grp_data']['vspan_half']

                    x = grid[:, :, 0]
                    y = grid[:, :, 1]
                    u = RF[:, :, 0]*-1
                    v = RF[:, :, 1]

                    if names=='FlpND':
                        alpha=0.75
                        z=2
                    elif names=='FlpD':
                        alpha=0.7
                        z=3
                    # z=2
                    # alpha=1

                    ax_main.quiver(x, y, u, v, width=0.005, scale=1.5, headwidth=4, headlength=6, headaxislength=5, color=color_dict[names][0],
                                alpha=alpha, minlength=1, zorder=z)
                    
                    if draw_contour:
                        ## draw contour for all responses
                        RF_amp_expanded = zoom(RF_amp, [expansion_factor, expansion_factor])
                        grid_new = expand_grid(13*expansion_factor, 9*expansion_factor)
                        x_new = grid_new[:, :, 0]
                        y_new = grid_new[:, :, 1]
                        print('here')
                        ax_main.contourf(x_new, y_new, RF_amp_expanded/np.amax(RF_amp_expanded), contour_levels, vmin=0, vmax=1, algorithm='serial',
                                        cmap=ListedColormap(create_color_gradient(color_dict[names][0], [1, 1, 1], 128)), alpha=0.2, antialiased=False)
                    else:
                        pass

                    if draw_contour_only_for_horizontal:
                        ## draw contour only for horizontal responses
                        RF_amp_expanded = zoom(RF_amp, [expansion_factor, expansion_factor])
                        grid_new = expand_grid(13*expansion_factor, 9*expansion_factor)
                        x_new = grid_new[:, :, 0]
                        y_new = grid_new[:, :, 1]
                        ax_main.contourf(x_new, y_new, u/np.amax(u), contour_levels, linewidths=0.4, vmin=0, vmax=1, 
                                        cmap=ListedColormap(create_color_gradient(color_dict[names][0], [1, 1, 1], 128)), alpha=0.1)
                    else:
                        pass
                    
                    if draw_grid_colors:
                        cmap = cm.LinearSegmentedColormap.from_list("mycmap", create_color_gradient(color_dict[names][0], [0.9, 0.9, 0.9], 128))
                        ax_main.pcolormesh(x, y, RF_amp, cmap=cmap, vmin=0, vmax=0.15)
                    else:
                        pass


                    x_axis = np.linspace(x[0, 0], x[0, -1], grid.shape[1])
                    ax_xDist.bar(x_axis, hor, width=x_axis[1]-x_axis[0], color=color_dict[names][0], alpha=0.5, linewidth=0.5, edgecolor='k')
                    y_axis = np.linspace(y[0, 6], y[-1, 6], grid.shape[0])
                    ax_yDist.barh(y_axis, vert, y_axis[1]-y_axis[0], color=color_dict[names][0], alpha=0.5, linewidth=0.5, edgecolor='k')

                ## for normal plotting
                # ax_xDist.spines['left'].set_visible(False)
                # ax_yDist.spines['bottom'].set_visible(False)
                # ax_xDist.spines['right'].set_visible(True)
                # ax_yDist.spines['top'].set_visible(True)
                ax_main.spines['left'].set_visible(False)
                ax_main.spines['bottom'].set_visible(False)
                ax_main.spines['right'].set_visible(False)
                ax_main.spines['top'].set_visible(False)
                ax_xDist.spines['left'].set_visible(False)
                ax_yDist.spines['bottom'].set_visible(False)
                ax_xDist.spines['right'].set_visible(True)
                ax_yDist.spines['top'].set_visible(True)
                # ax_yDist.spines['left'].set_visible(False)
                # ax_xDist.spines['bottom'].set_visible(False)
                ax_yDist.spines['left'].set_position(('data', 0))
                ax_xDist.spines['bottom'].set_position(('data', 0))
                # ax_xDist.spines['bottom'].set_bounds(low=x[0, -1], high=x[0, 0])
                # ax_yDist.spines['left'].set_bounds(low=y[-1, 6], high=y[0, 6])
                ##
                ax_xDist.yaxis.tick_right()
                ax_yDist.xaxis.tick_top()
                ax_main.yaxis.tick_left()
                ax_main.xaxis.tick_bottom()
                ##
                ax_xDist.set_ylim(top=response, bottom=min_response)
                ax_yDist.set_xlim(left=min_response, right=response)
                ax_xDist.set_yticks([response])
                ax_yDist.set_xticks([response])
                ax_xDist.spines['right'].set_bounds(low=min_response, high=response)
                ax_yDist.spines['top'].set_bounds(low=min_response, high=response)
                ##
                n=10
                ax_main.set_xlim(x[0, -1]+(x[0, -1]/n), x[0, 0]+(x[0, 0]/n))
                ax_xDist.set_xlim(x[0, -1]+(x[0, -1]/n), x[0, 0]+(x[0, 0]/n))
                # y_min = np.amin(y)
                # y_max = np.amax(y)
                # m=9
                # ax_main.set_ylim(y_min+(y_min/m), y_max+(y_max/m))
                # ax_yDist.set_ylim(y_min+(y_min/m), y_max+(y_max/m))
                ##
                ax_main.set_xticks([x[0, -1], 0, x[0, 0]])
                ax_main.set_yticks([y[-1, 6], 0, y[0, 6]])
                ax_xDist.set_xticks([])
                ax_yDist.set_yticks([])
                fig.savefig(os.path.join(Dir_output, filename+'_'+cellname+direction+'deformed_RFs.png'), dpi=600)
    return 1

## Fig4h,j
def compare_contra_and_ipsi_response(cellname_list=['HSN'], direction_list=['neg_x'], flies_to_plot_list=['FlpND', 'FlpD'], 
                                     Dir=r'C:\Users\rsatapat\Documents\Victoria\RepositoryForScanningRect\res'):
    total_df = pd.DataFrame()
    for cellname in cellname_list:
        for names in flies_to_plot_list:
            df = pd.DataFrame()
            data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname+'_' + names + '.mat'))['grp_data']['RFall']
            data_ipsi = np.linalg.norm(np.sum(data[:, :7, :2, :], axis=(0, 1)), axis=0)
            data_contra = np.linalg.norm(np.sum(data[:, 7:, :2, :], axis=(0, 1)), axis=0)
            data_ratio = np.divide(data_contra, data_ipsi)
            a=np.array([data_ipsi, data_contra, data_ratio]).T
            df = pd.concat([df, pd.DataFrame(a, columns=['ipsi', 'contra', 'ratio'])], axis=0)
            df = pd.melt(df, value_vars=['ipsi', 'contra', 'ratio'], var_name='side', value_name='response')
            df = pd.concat([df, pd.DataFrame(data={'strain':[names]*data.shape[-1]*3})], axis=1)
            total_df = pd.concat([total_df, df], axis=0, ignore_index=True)
  
        fig, ax = plt.subplots(figsize=(2, 3))
        palette = [shades_ND[5], shades_D[3], shades_DB331[2]]
        # ax.boxplot(list(total_df.loc[total_df['type_gen']=='VSshakB2']['mean-JP']), widths=0.3, positions=[3], showfliers=False, showbox=False, showcaps=False, showmeans=True,
        #               meanline=True, medianprops = dict(linewidth=0), meanprops = dict(linestyle='solid', linewidth=1, color=color_dict['ShakB2'][0]), whiskerprops=dict(linewidth=0), zorder=3)
        # ax.errorbar(x=[3], y=np.mean(list(total_df.loc[total_df['type_gen']=='VSshakB2']['mean-JP'])), yerr=sem(list(total_df.loc[total_df['type_gen']=='VSshakB2']['mean-JP'])),
        #              ecolor=color_dict['ShakB2'][0], elinewidth=1)
        # ratio(ipsi/contra)
        sns.pointplot(y=total_df.loc[total_df['side']=='ratio']['response'], x='strain', order=flies_to_plot_list, hue='strain',  hue_order=flies_to_plot_list, data=total_df, dodge=0.4, ax=ax, palette=palette, estimator='mean', errorbar='se', join=False,
                    markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.75, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.75, lw=1.5)
        sns.stripplot(y=total_df.loc[total_df['side']=='ratio']['response'], x='strain', order=flies_to_plot_list, hue='strain',  hue_order=flies_to_plot_list, data=total_df, dodge=0.5, ax=ax, size=10, palette=palette, alpha=0.55, edgecolor='k', linewidth=1, jitter=0.05)

        ax.set_ylim(top=1, bottom=-0.005)
        annotator = Annotator(ax, data=total_df, y=total_df.loc[total_df['side']=='ratio']['response'], x='strain', order=flies_to_plot_list, hue='strain',  hue_order=flies_to_plot_list,
                            pairs = [((flies_to_plot_list[0], flies_to_plot_list[0]), (flies_to_plot_list[i], flies_to_plot_list[i])) for i in range(1, len(flies_to_plot_list))])
                            # pairs=[(('FlpND', 'FlpND'), ('FlpD', 'FlpD')), (('FlpND', 'FlpND'), ('FlpNDxDB331', 'FlpNDxDB331'))])
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
        _, annotations = annotator.apply_and_annotate()

        ax.spines['bottom'].set_visible(False)
        ax.set_title('Total response ipsi contra')
        ax.set_ylabel('Response')
        ax.set_xlabel('')
        ax.tick_params(bottom=False)
        ax.legend_.remove()
        ax.spines['left'].set_bounds(low=0, high=1)
        ax.spines['left'].set_position(('outward', 10))
        ax.set_xlim(left=-0.5, right=2)

        fig.savefig(os.path.join(Dir_output, cellname +'contra_ipsi_ratio_response.png'))
        fig.savefig(os.path.join(Dir_output, cellname +'contra_ipsi_ratio_response.pdf'), transparent=True, dpi=600)

        with open(os.path.join(Dir_output, cellname +'contra_ipsi_response.txt'), 'w') as f:
            f.write('testname' + ',' + 'genotypes' + ',' + 'cut' + ',' + 'U-val' + ',' + 'p-val' + ',' + '#stars' + "\n")
            for i in range(len(annotations)):
                f.write(annotations[i].formatted_output.split(',')[0] +annotations[i].structs[0]['label']+ '_' +annotations[i].structs[1]['label'] +
                        ',' + annotations[i].formatted_output.split(',')[1].split(' ')[1].split(':')[1]+ ',' +
                        annotations[i].formatted_output.split(',')[1].split(' ')[2].split('=')[1] + ',' + str(len(annotations[i].text)) +"\n")
                
        sns.pointplot(y="response", x='side', order=['ipsi', 'contra'], hue='strain',  hue_order=flies_to_plot_list, data=total_df, dodge=0.4, ax=ax, palette=palette, estimator='mean', errorbar='se', join=False,
                markers="_", errwidth=1, scale=2)
        plt.setp(ax.lines, zorder=100, alpha=0.75, lw=1.5)
        plt.setp(ax.collections, zorder=100, alpha=0.75, lw=1.5)
        sns.stripplot(y="response", x='side', order=['ipsi', 'contra'], hue='strain',  hue_order=flies_to_plot_list, data=total_df, dodge=0.5, ax=ax, size=10, palette=palette, alpha=0.55, edgecolor='k', linewidth=1, jitter=0.05)

        ax.set_ylim(top=7, bottom=-0.005)
        annotator = Annotator(ax, data=total_df, y="response", x='side', order=['ipsi', 'contra'], hue='strain', hue_order=flies_to_plot_list, #pairs=[('FlpND', 'FlpD')])
                            pairs = [(('ipsi', flies_to_plot_list[0]), ('ipsi', flies_to_plot_list[i])) for i in range(1, len(flies_to_plot_list))] + 
                            [(('contra', flies_to_plot_list[0]), ('contra', flies_to_plot_list[i])) for i in range(1, len(flies_to_plot_list))])

        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
        _, annotations = annotator.apply_and_annotate()

        ax.spines['bottom'].set_visible(False)
        ax.set_title('Total response ipsi contra')
        ax.set_ylabel('Response')
        ax.set_xlabel('')
        ax.tick_params(bottom=False)
        ax.legend_.remove()
        ax.spines['left'].set_bounds(low=0, high=7)
        ax.spines['left'].set_position(('outward', 10))
        ax.set_xlim(left=-0.5, right=2)

        fig.savefig(os.path.join(Dir_output, cellname +'contra_ipsi_response.png'))
        fig.savefig(os.path.join(Dir_output, cellname +'contra_ipsi_response.svg'), transparent=True, dpi=600)

        with open(os.path.join(Dir_output, cellname +'contra_ipsi_response.txt'), 'w') as f:
            for i in range(len(annotations)):
                f.write(annotations[i].formatted_output+"\n")
    return 1

## Fig4e,f
def plot_rawRF_difference_as_dots(flies_to_plot_list=['FlpND', 'FlpD'], cellname_list=['HSN'], direction='neg_x'):
    RF_data = []
    RFx_raw = []
    RFy_raw = []
    for cellname in cellname_list:
        for names in flies_to_plot_list:
            data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname + '_' + names + '.mat'))['grp_data']
            RF_data.append(data)
            RFx_raw.append(data['RF_x_pos_raw'] - data['RF_x_neg_raw'])
            RFy_raw.append(data['RF_y_pos_raw'] - data['RF_y_neg_raw'])
            print(np.amax(RFx_raw[-1]), np.amin(RFx_raw[-1]))

    if 'pos' in direction:
        RF_diff = np.subtract(RF_data[0]['RF'+'_'+direction+'_raw']*-1, RF_data[1]['RF'+'_'+direction+'_raw']*-1)
    elif 'neg' in direction:
        RF_diff = np.subtract(RF_data[0]['RF' + '_' + direction + '_raw'] * -1, RF_data[1]['RF' + '_' + direction + '_raw'] * -1) * -1
    else:
        RF_diff = np.subtract(RFx_raw[0]*-1, RFx_raw[1]*-1)
        # RFy_diff = np.subtract(RFy_raw[0]*-1, RFy_raw[1]*-1)
        # RFx_sum = np.add(RFx_raw[0]*-1, RFx_raw[1]*-1)
        # RFy_sum = np.add(RFy_raw[0]*-1, RFy_raw[1]*-1)
        # RF_diff = np.divide(np.linalg.norm(np.stack((RFx_diff, RFy_diff), axis=2), axis=2),
        #                     np.linalg.norm(np.stack((RFx_sum, RFy_sum), axis=2), axis=2))
    grid = RF_data[0]['deformed_grid']

    x = grid[:, :, 0]
    y = grid[:, :, 1]
    cmap = cm.LinearSegmentedColormap.from_list("mycmap", create_color_gradient(color_dict['FlpD'][0], [0.99, 0.99, 0.99], 128)[::-1] +
                                                create_color_gradient(color_dict['FlpND'][0], [0.99, 0.99, 0.99], 128))
    print(np.amax(RF_diff), np.min(RF_diff), np.amin(RF_diff))
    # normalize = matplotlib.colors.Normalize(vmin=-np.amax(RF_diff)/2, vmax=np.amax(RF_diff)/2)
    normalize = matplotlib.colors.Normalize(vmin=-0.014/2, vmax=0.014/2)
    colors = cmap(normalize(np.array(RF_diff).flatten()))
    fig, ax = plt.subplots()
    print(x.shape, y.shape, RF_diff.shape)
    #np.square(RF_diff)*2500000
    ax.scatter(x, y, s=np.square(RF_diff)*500000, c=colors, alpha=0.7, linewidths=0)

    p1 = ax.scatter(-0.8, 0.65, s=np.square(0.015)*500000, c='k', alpha=0.7, linewidths=0, label='1.5')
    p2 = ax.scatter(-0.8, 0.63, s=np.square(0.010)*500000, c='k', alpha=0.7, linewidths=0, label='1.0')
    p3 = ax.scatter(-0.8, 0.61, s=np.square(0.005)*500000, c='k', alpha=0.7, linewidths=0, label='0.5')
    ax.legend(loc='upper right', ncol=1, fontsize=8, frameon=False)
    p1.remove()
    p2.remove()
    p3.remove()
    n=10
    ax.set_xlim(x[0, -1]+(x[0, -1]/n), x[0, 0]+(x[0, 0]/n))
    ax.set_xticks([x[0, -1], 0, x[0, 0]])
    ax.set_yticks([y[-1, 6], 0, y[0, 6]])
    ax.set_title(cellname_list[0])
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    fig.savefig(os.path.join(Dir_output, flies_to_plot_list[0]+'_'+flies_to_plot_list[1]+cellname_list[0]+'_raw_diff.png'), dpi=300)
    fig.savefig(os.path.join(Dir_output, flies_to_plot_list[0]+'_'+flies_to_plot_list[1]+cellname_list[0]+'_raw_diff.pdf'), dpi=300)

    # flies_to_plot_list = [x.replace('_spike', '') if 'spike' in x else x for x in flies_to_plot_list]
    # fig, ax = plt.subplots()
    # ax.bar(np.arange(0, 13), np.mean(RF_data[0]['RF'+'_'+direction+'_raw']*-1, axis=0), color=color_dict[flies_to_plot_list[0]][0], alpha=0.5)
    # ax.bar(np.arange(0, 13), np.mean(RF_data[1]['RF'+'_'+direction+'_raw']*-1, axis=0), color=color_dict[flies_to_plot_list[1]][0], alpha=0.5)
    # n=10
    # ax.set_xlim(left=13, right=-1)
    # fig.savefig(os.path.join(Dir_output, flies_to_plot_list[0]+'_'+flies_to_plot_list[1]+'_ND_comparison.png'), dpi=300)
    
    return 1

## supplementary fig 11f,g
def plot_normalized_azimuth_and_elevation_response(cellname = 'HSE', flies_to_plot = ['FlpND', 'FlpD'], filename = 'pos_x', direction = '',
                                                   Dir=r'C:\Users\rsatapat\Documents\Victoria\RepositoryForScanningRect\res'):
    fig, ax = plt.subplots(figsize=(4, 3))
    dict = {}
    total_azimuth_data  = []
    avg_contra_ipsi_comparison = []
    for names in flies_to_plot:
        filename += names
        data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname+'_' + names + '.mat'))

        azimuth_resp_normed = data['grp_data']['horizontal_profiles_normed1']
        elevation_resp_normed = data['grp_data']['vertical_profiles_normed1']
        total_azimuth_data.append(azimuth_resp_normed)

        ## for each cell
        total_area_eachcell = np.apply_along_axis(lambda x: np.trapz(x, dx=1), axis=1, arr=azimuth_resp_normed)
        ipsi_area_eachcell = np.apply_along_axis(lambda x: np.trapz(x, dx=1), axis=1, arr=azimuth_resp_normed[:, 6:])
        contra_area_eachcell = np.apply_along_axis(lambda x: np.trapz(x, dx=1), axis=1, arr=azimuth_resp_normed[:, :7])
        dict[names] = {'total': total_area_eachcell, 'ipsi': ipsi_area_eachcell, 'contra': contra_area_eachcell}

        # azimuth_resp_normed_interp = []
        # for i in range(azimuth_resp_normed.shape[0]):
        #     interp = interp1d(np.arange(0, azimuth_resp_normed.shape[1]), azimuth_resp_normed[i], kind='cubic')
        #     azimuth_resp_normed_interp.append(interp(np.linspace(0, 12, 30)))
        # azimuth_resp_normed = np.array(azimuth_resp_normed_interp)
        azimuth_mean_response = np.mean(azimuth_resp_normed, axis=0)
        azimuth_sem_response = sem(azimuth_resp_normed, axis=0)
        elevation_mean_response = np.mean(elevation_resp_normed, axis=0) * 5
        elevation_sem_response = sem(elevation_resp_normed, axis=0) * 5

        ## for average
        total_area = np.trapz(azimuth_mean_response, dx=1)
        ipsi_area = np.trapz(azimuth_mean_response[6:], dx=1)
        contra_area = np.trapz(azimuth_mean_response[:7], dx=1)
        print(total_area/total_area, ipsi_area/total_area, contra_area/total_area)
        avg_contra_ipsi_comparison.append([ipsi_area/total_area, contra_area/total_area])

        # plot each neuron
        # for i in range(azimuth_resp_normed.shape[0]):
        #     ax.plot(np.arange(0, 13), azimuth_resp_normed[i], color=color_dict[names][0], alpha=0.2, label=names, lw=0.3)
        ax.plot(np.arange(0, 13), azimuth_mean_response, color=color_dict[names][0], alpha=0.8, label=names, lw=2)
        ax.fill_between(np.arange(0, 13), y1=azimuth_mean_response-azimuth_sem_response, y2=azimuth_mean_response+azimuth_sem_response,
                                    color=color_dict[names][0], alpha=0.1)

        # for plotting both azimuth and elevation in the same plot
        # ax2 = ax.twinx()
        # ax2.plot(elevation_mean_response, np.arange(0, 9), color=color_dict[names][0], alpha=0.8, label=names)
        # ax2.fill_betweenx(np.arange(0, 9), x1=elevation_mean_response-elevation_sem_response, x2=elevation_mean_response+elevation_sem_response,
        #                             color=color_dict[names][0], alpha=0.1)

    np.save(os.path.join(Dir, 'RF_azimuth_' + cellname + '_area.npy'), dict, allow_pickle=True)
    area_data = np.load(os.path.join(Dir, 'RF_azimuth_' + cellname + '_area.npy'), allow_pickle=True)[()]
    #
    textfile = os.path.join(Dir_output, 'RF_azimuth_'+cellname+'_significance_test.txt')
    testname = 'Mann-WhitneyU'
    with open(textfile, 'w') as f:
        f.write('testname' + '       ' + 'U-val' + '   ' + 'p-val  ' + 'stars'+"\n")
        for i in range(total_azimuth_data[0].shape[1]):
            test = mannwhitneyu(total_azimuth_data[0][:, i], total_azimuth_data[1][:, i])
            stars = significance([test.pvalue])[0]
            f.write(testname + '   ' + str(test.statistic) + '   ' + str(round(test.pvalue, 4))+'   '+ str(stars)+"\n")
    f.close()
    #
    textfile = os.path.join(Dir_output, 'RF_azimuth_' + cellname + '_area_significance_test.txt')
    testname = 'Mann-WhitneyU'
    with open(textfile, 'w') as f:
        f.write('testname' + '       ' +'type   '+ 'U-val' + '   ' + 'p-val   ' + 'stars' "\n")
        test = mannwhitneyu(area_data['FlpND']['ipsi']/area_data['FlpND']['total'], area_data['FlpD']['ipsi']/area_data['FlpD']['total'])
        stars = significance([test.pvalue])[0]
        f.write(testname+'   '+ 'ipsi   '+str(test.statistic)+'   '+str(round(test.pvalue, 4))+'   '+ str(stars)+"\n")
        test2 = mannwhitneyu(area_data['FlpND']['contra']/area_data['FlpND']['total'], area_data['FlpD']['contra']/area_data['FlpD']['total'])
        stars = significance([test2.pvalue])[0]
        f.write(testname+'   '+ 'contra '+str(test2.statistic)+'   '+str(round(test2.pvalue, 4))+'   '+ str(stars)+"\n")
        f.write('FlpND   '+'ipsi :'+ str(round(avg_contra_ipsi_comparison[0][0], 4)*100)+'  contra :'+
                str(round(avg_contra_ipsi_comparison[0][1], 4)*100)+'\n')
        f.write('FlpD    '+'ipsi :'+ str(round(avg_contra_ipsi_comparison[1][0], 4)*100)+'  contra :'+
                str(round(avg_contra_ipsi_comparison[1][1], 4)*100)+'\n')
    f.close()

    ax.set_ylim(top=1, bottom=0)
    ax.set_xlim(12.2, -0.2) ## reversed to reverse the x-axis
    # ax.set_xlim(10.2, -0.2) ## reversed to reverse the x-axis
    ax.spines['left'].set_bounds(low=0, high=1)
    ax.spines['bottom'].set_bounds(low=0, high=12)
    # ax.spines['bottom'].set_bounds(low=0, high=10)
    ax.set_yticks([0, 1])
    ax.set_xticks([0, 6, 12])
    # ax.set_xticks([0, 5, 10])
    ax.set_xticklabels(['70', '0', '-70'])
    ax.axvline(x=6, ls='--', lw=1, c='grey')
    # ax.axvline(x=5, ls='--', lw=1, c='grey')
    
    # ax2.xaxis.tick_top()
    # ax2.yaxis.tick_right()
    # ax2.spines['right'].set_visible(True)
    # # ax2.spines['top'].set_visible(True)
    # ax2.spines['left'].set_visible(True)
    # ax2.spines['bottom'].set_visible(True)
    # ax2.set_ylim(-0.2, 8.2)
    # ax.legend()
    fig.savefig(os.path.join(Dir_output, filename+'_'+cellname+direction+'_az_resp.png'), dpi=300)
    fig.savefig(os.path.join(Dir_output, filename+'_'+cellname+direction+'_az_resp.pdf'), dpi=600)
    
    return 1

## Supplementary Fig 11h,i
def plot_contour_areas(cellname='HSE', flies_to_plot=['FlpND', 'FlpD'], Dir=r'C:\Users\rsatapat\Documents\Victoria\RepositoryForScanningRect\res'):    
    fig, ax = plt.subplots(figsize=(5, 3))
    total_cnt_areas = []
    for i in range(len(flies_to_plot)):
        data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname+'_' + 'FlpD' + '.mat'))
        cnt_levels = data['grp_data']['cntr_levels']
        cnt_areas = data['grp_data']['cntr_areas']
        total_cnt_areas.append(cnt_areas)
        for i in range(cnt_areas.shape[0]):
            ax.plot(cnt_levels, cnt_areas[i], c=color_dict[flies_to_plot][0], lw=0.5, alpha=0.25)

        ax.plot(cnt_levels, np.mean(cnt_areas, axis=0), c=color_dict[flies_to_plot][0], lw=2)
    
    ax.spines['left'].set_bounds(low=1, high=0)
    ax.spines['bottom'].set_bounds(low=cnt_levels[0], high=cnt_levels[-1])
    ax.set_xlim(left=0.08)
    ax.set_xticks(cnt_levels)
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.set_xticklabels([str(int(x*100)) for x in cnt_levels])
    print(cnt_levels)
    ax.grid(visible=True, alpha=0.7, lw=0.5, ls='dashed')
    textfile = os.path.join(Dir_output, 'cntr_area_'+cellname+'_significance_test.txt')
    testname = 'Mann-WhitneyU'
    with open(textfile, 'w') as f:
        f.write('testname' + '       ' + 'value' + ' ' + 'U-val' + '   ' + 'p-val' + "\n")
        for j in range(1, len(total_cnt_areas)):
            for i in range(total_cnt_areas[0].shape[1]):
                test = mannwhitneyu(total_cnt_areas[0][:, i], total_cnt_areas[i][:, i])
                f.write(testname+'   '+str(int(cnt_levels[i]*100))+'   '+str(test.statistic)+'   '+str(round(test.pvalue, 4))+"\n")
    f.close()
    fig.savefig(os.path.join(Dir_output, 'cntr_area_'+cellname+'.png'), dpi=300)
    fig.savefig(os.path.join(Dir_output, 'cntr_area_'+cellname+'.pdf'), dpi=600)

    return 1
    
## Fig 4g,i
def plot_contour_areas_at_cutoff(cellname='HSE', flies_to_plot=['FlpND', 'FlpD'], Dir=r'C:\Users\rsatapat\Documents\Victoria\RepositoryForScanningRect\res'):
    # bar plot comparing size at 50% cutoff, significance test
    fig, ax = plt.subplots(figsize=(2, 3))

    total_cnt_areas = []
    for i in range(len(flies_to_plot)):
        data = mat73.loadmat(os.path.join(Dir, 'scanning_rect_grp_' + cellname+'_' + 'FlpD' + '.mat'))
        cnt_levels = data['grp_data']['cntr_levels']
        cnt_areas = data['grp_data']['cntr_areas']
        total_cnt_areas.append(cnt_areas)

    cutoff = 50
    cnt_levels = [str(int(x * 100)) for x in cnt_levels]
    cnt_areas_total = np.concatenate(total_cnt_areas, axis=0)
    df = pd.DataFrame(data=cnt_areas_total, columns=cnt_levels)
    print([[flies_to_plot[0]]*total_cnt_areas[0].shape[0] + [flies_to_plot[1]]*total_cnt_areas[1].shape[0]])
    df['strain'] = [flies_to_plot[0]]*total_cnt_areas[0].shape[0] + [flies_to_plot[1]]*total_cnt_areas[1].shape[0]

    data = df
    data_melted = pd.melt(data, id_vars='strain', var_name='cutoffs', value_name='cntr_area')

    palette = [shades_ND[5], shades_D[3]]
    # plot for multiple cutoffs
    sns.pointplot(y='cntr_area', x='cutoffs', order=['30', '60'], hue='strain', hue_order=flies_to_plot, data=data_melted, dodge=0.4, ax=ax, palette=palette, estimator='mean', errorbar='se', join=False,
                  markers="_", errwidth=1, scale=2)
    plt.setp(ax.lines, zorder=100, alpha=0.95, lw=1.5)
    plt.setp(ax.collections, zorder=100, alpha=0.95, lw=1.5)
    sns.stripplot(y='cntr_area', x='cutoffs', order=['30', '60'], hue='strain', hue_order=flies_to_plot, data=data_melted, dodge=True, ax=ax, size=7, palette=palette, alpha=0.65, edgecolor='k',
                  linewidth=0.4, jitter=0.25)
    
    # ax.set_ylim(top=0.7, bottom=-0.005)
    annotator = Annotator(ax, data=data_melted, y='cntr_area', x='cutoffs', order=['30', '60'], hue='strain', hue_order=flies_to_plot,
                          pairs = [(('30', flies_to_plot[0]), ('30', flies_to_plot[i])) for i in range(1, len(flies_to_plot))] +
                                    [(('60', flies_to_plot[0]), ('60', flies_to_plot[i])) for i in range(1, len(flies_to_plot))])
    annotator.configure(test='Mann-Whitney-gt', text_format='star', loc='inside')
    _, annotations = annotator.apply_and_annotate()
    
    ax.set_xlim(left=-0.5, right=1.5)
    ax.set_yticks(np.arange(0, 1.05, 0.2))
    ax.set_yticklabels([str(int(x * 100)) for x in np.arange(0, 1.05, 0.2)])
    ax.spines['bottom'].set_visible(False)
    ax.set_title(cellname + ' RF size at 30% and 60% cutoff')
    ax.set_ylabel('RF area (% of total)')
    ax.set_xlabel('')
    ax.tick_params(bottom=False)
    ax.legend_.remove()
    ax.spines['left'].set_bounds(low=0, high=1)
    
    fig.savefig(os.path.join(Dir_output, 'RF_areas'+cellname+'.png'))
    fig.savefig(os.path.join(Dir_output, 'RF_areas'+cellname+'.pdf'), transparent=True, dpi=600)
    
    with open(os.path.join(Dir_output, 'RF_areas'+cellname+'.txt'), 'w') as f:
        for i in range(len(annotations)):
            f.write(annotations[i].formatted_output+"\n")
    f.close()
    
    # plot for just one cutoff
    sns.pointplot(y=str(cutoff), x='strain', hue='strain', hue_order=flies_to_plot, data=data, dodge=0.4, ax=ax, palette=palette, estimator='mean', errorbar='se', join=False,
                  markers="_", errwidth=1, scale=2)
    plt.setp(ax.lines, zorder=100, alpha=0.75, lw=1.5)
    plt.setp(ax.collections, zorder=100, alpha=0.75, lw=1.5)
    sns.stripplot(y=str(cutoff), x='strain', hue='strain', hue_order=flies_to_plot, data=data, dodge=True, ax=ax, size=10, palette=palette, alpha=0.55, edgecolor='k', linewidth=0.4, jitter=0.05)
    
    ax.set_ylim(top=0.7, bottom=-0.005)
    annotator = Annotator(ax, data=data, y=str(cutoff), x='strain', order=flies_to_plot, hue='strain', hue_order=flies_to_plot,
                          pairs=[((flies_to_plot[0], flies_to_plot[0]), (flies_to_plot[i], flies_to_plot[i])) for i in range(1, len(flies_to_plot))])
                        #   pairs=[(('FlpND', 'FlpND'), ('FlpD', 'FlpD'))])
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
    _, annotations = annotator.apply_and_annotate()
    
    ax.set_xlim(left=-0.5, right=1.5)
    ax.set_yticks(np.arange(0, 0.65, 0.2))
    ax.set_yticklabels([str(int(x * 100)) for x in np.arange(0, 0.65, 0.2)])
    ax.spines['bottom'].set_visible(False)
    ax.set_title('RF size at '+str(cutoff)+'% cutoff')
    ax.set_ylabel('RF area (% of total)')
    ax.set_xlabel('')
    ax.tick_params(bottom=False)
    ax.legend_.remove()
    ax.spines['left'].set_bounds(low=0, high=0.6)
    
    fig.savefig(os.path.join(Dir_output, 'RF_areas'+cellname+str(cutoff)+'.png'))
    fig.savefig(os.path.join(Dir_output, 'RF_areas'+cellname+str(cutoff)+'.svg'), transparent=True, dpi=600)
    
    with open(os.path.join(Dir_output, 'RF_areas'+cellname+str(cutoff)+'.txt'), 'w') as f:
        for i in range(len(annotations)):
            f.write(annotations[i].formatted_output+"\n")
    f.close()
    return 1