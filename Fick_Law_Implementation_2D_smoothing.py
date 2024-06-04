# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 17:55:40 2021
@author: Ashlyn
"""
from matplotlib.colors import LinearSegmentedColormap
import datetime
import numpy as np
import csv
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from numpy import trapz

'---Image settings-----'
size = 20
plt.rc('axes', titlesize=size)     # fontsize of the axes title
plt.rc('axes', labelsize=size)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=size)    # fontsize of the tick labels
plt.rc('legend', fontsize=size)      # legend fontsize
plt.rc('figure', titlesize=5)      # fontsize of the figure title
mpl.rcParams['figure.figsize'] = 4, 5
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10

'define RMSE'
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

'define colormap'
def cmap_generation(length):
    alpha = np.linspace(1, 0.4, length)
    cmap = LinearSegmentedColormap.from_list(
        name='camp', colors=['xkcd:robin egg blue', 'magenta'], N=length)
    colors = []
    for ind in range(cmap.N):
        c = []
        aa = alpha[ind]
        for x in cmap(ind)[:3]:
            c.append(x*aa)
        colors.append(tuple(c))
    my_cmap = mpl.colors.ListedColormap(colors)
    return my_cmap

save_file_path_1 = 'savefig/02_Ficks_Law_Implementation/'
save_file_path_2 = 'savefig/03_Ficks_Law_Concentration/'
# save_file_path_3 = 'savefig/02_Ficks_Law_Implementation/003'
# save_fig_1 = True  # DNA vs z-direction
save_fig_2 = True  # RA vs DF & Amomalous map
save_fig_3 = True   # RA vs DF & Amomalous map + Whisker Outliers

Nz = 100
Ns = 2
# file_path = 'Crystals_smooth_csv/'
file_path = 'new_Crystal_smooth_csv/'

'Crystal 1'
name01 = 'Crystal1_fold'
Lp1 = 27.2*10**-4

'Crystal 2'
name02 = 'Crystal2_fold'
Lp2 = 20.4*10**-4

'Crystal 3'
name03 = 'Crystal3_fold'
Lp3 = 61.88*10**-4

all_Crystal_num = [1.0, 2.0, 3.0]
all_name123 = [name01, name02, name03]
all_Lp123 = [Lp1, Lp2, Lp3]
all_xtal123 = ['C1', 'C2', 'C3']

all_Crystal_num = [1.0]
all_name123 = [name01]
all_Lp123 = [Lp1]
all_xtal123 = ['C1']

# all_Crystal_num = [2.0]
# all_name123 = [name02]
# all_Lp123 = [Lp2]
# all_xtal123 = ['C2']

# all_Crystal_num = [3.0]
# all_name123 = [name03]
# all_Lp123 = [Lp3]
# all_xtal123 = ['C3']

_ = 0
for xtal_ii, namei, xtal, Lp in zip(all_Crystal_num, all_name123, all_xtal123, all_Lp123):
    df1l = pd.DataFrame()
    df2l = pd.DataFrame()
    df3l = pd.DataFrame()
    df4l = pd.DataFrame()

    df1r = pd.DataFrame()
    df2r = pd.DataFrame()
    df3r = pd.DataFrame()
    df4r = pd.DataFrame()
    csv_index = -1
    _ += 1
    file_path_2 = 'new_Whisker_plot_csv/' 
    
    filename1_l = file_path_2+'all_crystal_20fold_smooth_c{:d}_Dl.csv'.format(_)
    filename2_l = file_path_2+'all_crystal_20fold_smooth_c{:d}_Timel.csv'.format(_)
    filename3_l = file_path_2+'all_crystal_20fold_smooth_c{:d}_Concl.csv'.format(_)
    filename4_l = file_path_2+'all_crystal_20fold_smooth_c{:d}_Depthl.csv'.format(_)

    filename1_r = file_path_2+'all_crystal_20fold_smooth_c{:d}_Dr.csv'.format(_)
    filename2_r = file_path_2+'all_crystal_20fold_smooth_c{:d}_Timer.csv'.format(_)
    filename3_r = file_path_2+'all_crystal_20fold_smooth_c{:d}_Concr.csv'.format(_)
    filename4_r = file_path_2+'all_crystal_20fold_smooth_c{:d}_Depthr.csv'.format(_)
    with open(filename1_l, 'w') as csvfile:
        pass
    with open(filename2_l, 'w') as csvfile:
        pass
    with open(filename3_l, 'w') as csvfile:
        pass
    with open(filename4_l, 'w') as csvfile:
        pass
    with open(filename1_r, 'w') as csvfile:
        pass
    with open(filename2_r, 'w') as csvfile:
        pass
    with open(filename3_r, 'w') as csvfile:
        pass
    with open(filename4_r, 'w') as csvfile:
        pass

    for xtal_fold_i in range(10):
        name1 = namei + str(xtal_fold_i)
        table_filename_l = 'Transport_tables/' + name1 + '_TransportTable_l.csv'
        with open(table_filename_l, 'w') as csvfile:
            csvfile.write(
                'dI/dz, Integral_trapz, Conc.(mol/cm^3), D(cm^2/sec), Time(sec), Position(element), depth(\u03BCm), crystal\n')

        table_filename_r = 'Transport_tables/' + name1 + '_TransportTable_r.csv'
        with open(table_filename_r, 'w') as csvfile:
            csvfile.write(
                'dI/dz, Integral_trapz, Conc.(mol/cm^3), D(cm^2/sec), Time(sec), Position(element), depth(\u03BCm), crystal\n')

        Crystal_num = [xtal_ii]
        name123 = [name1]
        Lp123 = [Lp]
        xtal123 = [xtal]

        if xtal == 'C1':
            t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
                      1087, 1207, 1327, 1447, 1567, 1687, 1807]
            # if save_fig_1 == True:
            #     mid_cutoff = 0
            # else:
            #     mid_cutoff = 0#10

        elif xtal == 'C2':
            t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
                      1087, 1207, 1327, 1447, 1567, 1687, 1807]
            # if save_fig_1 == True:
            #     mid_cutoff = 0
            # else:
            #     mid_cutoff = 0#5
        elif xtal == 'C3':
            t_eval = [17, 137, 257, 377, 497, 617, 737, 857, 977,
                      1097, 1217, 1337, 1457, 1577, 1697, 1817]
            # if save_fig_1 == True:
            #     mid_cutoff = 0
            # else:
            #     mid_cutoff = 0#5
        "Get image intensity from image detect"

        table_data_l = np.zeros([1, 8])
        table_data_r = np.zeros([1, 8])
        today = datetime.date.today()

        table_filename123_l = pd.DataFrame()
        table_filename123_r = pd.DataFrame()

        for xtal_i, name, xtal, Lp in zip(Crystal_num, name123, xtal123, Lp123):
            crystal = xtal_i
            print('File name=', name)
            Images_filename = file_path + name + '_images_smooth_data.csv'
            # Images_filename = name + '_images_FRAP_smoothestdata.csv'
            Images_data = pd.read_csv(Images_filename, delimiter=',')
            Nz = np.shape(Images_data)[1]

            for i in np.arange(len(t_eval)-1):

                x_array = np.linspace(0, Lp*10**4, Nz, endpoint=True)
                Integral_trapz_l_1 = 0
                dI_dz_l_1 = 0
                Integral_trapz_r_1 = 0
                dI_dz_r_1 = 0

                for j in np.arange(int(Nz/2 - 1)):
                    a = i+1  # time point: 1 - 14
                    b = j+1  # z-position: 1 - 49 (xtal edge to middle)

                    Integral_trapz_t2_l = trapz(
                        y=Images_data.iloc[a, b:int(Nz/2)], dx=(Lp/Nz))
                    Integral_trapz_t1_l = trapz(
                        y=Images_data.iloc[a-1, b:int(Nz/2)], dx=(Lp/Nz))
                    Integral_trapz_l = (Integral_trapz_t2_l-Integral_trapz_t1_l) / \
                        (t_eval[a]-t_eval[a-1])

                    Integral_trapz_t2_r = trapz(
                        y=Images_data.iloc[a, int(Nz/2):Nz-b], dx=(Lp/Nz))
                    Integral_trapz_t1_r = trapz(
                        y=Images_data.iloc[a-1, int(Nz/2):Nz-b], dx=(Lp/Nz))
                    Integral_trapz_r = (Integral_trapz_t2_r-Integral_trapz_t1_r) / \
                        (t_eval[a]-t_eval[a-1])

                    dI_dz_a_l = (
                        Images_data.iloc[a, b-1]-Images_data.iloc[a, b])/(Lp/Nz)
                    dI_dz_a_1_l = (
                        Images_data.iloc[a-1, b-1]-Images_data.iloc[a-1, b])/(Lp/Nz)
                    dI_dz_l = (dI_dz_a_l + dI_dz_a_1_l)/2

                    dI_dz_a_r = (
                        Images_data.iloc[a, (Nz-1)-b]-Images_data.iloc[a, (Nz-1)-(b-1)])/(Lp/Nz)
                    dI_dz_a_1_r = (
                        Images_data.iloc[a-1, (Nz-1)-b]-Images_data.iloc[a-1, (Nz-1)-(b-1)])/(Lp/Nz)
                    dI_dz_r = -(dI_dz_a_r + dI_dz_a_1_r)/2

                    Il = Images_data.iloc[a, b]
                    Conc_l = Il*845351.18
                    time_pt = t_eval[a-1]
                    depth_l = b*Lp/Nz*1e4  # um
                    position_l = b

                    Ir = Images_data.iloc[a, (Nz-1)-b]
                    Conc_r = Ir*845351.18
                    time_pt = t_eval[a-1]
                    depth_r = b*Lp/Nz*1e4  # um
                    position_r = Nz-b

                    Dl = (Integral_trapz_l-Integral_trapz_l_1)/(dI_dz_l-dI_dz_l_1)
                    Dr = (Integral_trapz_r-Integral_trapz_r_1)/(dI_dz_r-dI_dz_r_1)
                    Integral_trapz_l_1 = Integral_trapz_l
                    dI_dz_l_1 = dI_dz_l
                    Integral_trapz_r_1 = Integral_trapz_r
                    dI_dz_r_1 = dI_dz_r

                    table_data_l[0, 0] = dI_dz_l
                    table_data_l[0, 1] = Integral_trapz_l
                    table_data_l[0, 2] = Dl
                    table_data_l[0, 3] = Conc_l
                    table_data_l[0, 4] = time_pt
                    table_data_l[0, 5] = position_l
                    table_data_l[0, 6] = depth_l
                    table_data_l[0, 7] = crystal

                    csvfile = open(table_filename_l, 'a')
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerows(table_data_l)
                    csvfile.flush()

                    table_data_r[0, 0] = dI_dz_r
                    table_data_r[0, 1] = Integral_trapz_r
                    table_data_r[0, 2] = Dr
                    table_data_r[0, 3] = Conc_r
                    table_data_r[0, 4] = time_pt
                    table_data_r[0, 5] = position_r
                    table_data_r[0, 6] = depth_r
                    table_data_r[0, 7] = crystal

                    csvfile = open(table_filename_r, 'a')
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerows(table_data_r)
                    csvfile.flush()

            "Get data from table"
            tables_data_l = pd.read_csv(table_filename_l, delimiter=',')
            tables_data_l.columns = ['dI_dz', 'Integral_trapz', 'D',
                                     'Conc', 'time_pt', 'position', 'depth', 'crystal']
            tables_data_r = pd.read_csv(table_filename_r, delimiter=',')
            tables_data_r.columns = ['dI_dz', 'Integral_trapz', 'D',
                                     'Conc', 'time_pt', 'position', 'depth', 'crystal']

            table_filename123_l = pd.concat(
                [table_filename123_l, tables_data_l])
            table_filename123_r = pd.concat(
                [table_filename123_r, tables_data_r])

        tables_data_l = table_filename123_l
        tables_data_r = table_filename123_r
        crystal_names = {1.0: 'Crystal 1',
                         2.0: 'Crystal 2', 3.0: 'Crystal 3'}

        tables_data_l['crystal'] = tables_data_l['crystal'].map(crystal_names)
        tables_data_l.reset_index(inplace=True)

        tables_data_r['crystal'] = tables_data_r['crystal'].map(crystal_names)
        tables_data_r.reset_index(inplace=True)

        for i, c in zip(range(len(Crystal_num)), Crystal_num):
            Crystal_num[i] = crystal_names[c]

        'Normalize Data'
        poly_degree = 1
        fit_intercept = True
        show_fitting_result = 'No'
        show_fitting_result = 'Yes'

        'train_test_split data'
        target_name = 'D'
        numerical_columns = ['Conc', 'crystal']

        if xtal_fold_i in [1, 2, 3, 4, 5, 6, 7, 8]:
            csv_index += 1
            
            df1l.insert(csv_index, 'D', tables_data_l['D'], True)
            df2l.insert(csv_index, 'time', tables_data_l['time_pt'], True)
            df3l.insert(csv_index, 'Conc', tables_data_l['Conc'], True)
            df4l.insert(csv_index, 'Dpeth', tables_data_l['depth'], True)

            df1r.insert(csv_index, 'D', tables_data_r['D'], True)
            df2r.insert(csv_index, 'time', tables_data_r['time_pt'], True)
            df3r.insert(csv_index, 'Conc', tables_data_r['Conc'], True)
            df4r.insert(csv_index, 'Dpeth', tables_data_r['depth'], True)
        
        if xtal_fold_i == 4:
            'Left: RA vs DF: Conc'
            Title = 'Left RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            plot_data_l = tables_data_l
            plot_data_r = tables_data_r
            
            if xtal == 'C1':
                low_cut = 0.02
                high_cut = 0.035
                plot_data_conc_l = tables_data_l
                plot_data_conc_r = tables_data_r
            elif xtal == 'C2':
                low_cut = 0.007
                high_cut = 0.015
                plot_data_conc_l = tables_data_l
                plot_data_conc_r = tables_data_r
            elif xtal == 'C3':
                low_cut = 0.025
                high_cut = 0.05
                plot_data_conc_l = tables_data_l
                plot_data_conc_r = tables_data_r
            
            plot_data_l_0 = plot_data_conc_l.loc[plot_data_conc_l['Conc'] <= low_cut]
            plot_data_r_0 = plot_data_conc_r.loc[plot_data_conc_r['Conc'] <= low_cut]
            
            plot_data_l_1 = plot_data_conc_l.loc[(plot_data_conc_l['Conc'] > low_cut) & (plot_data_conc_l['Conc'] < high_cut)]
            plot_data_r_1 = plot_data_conc_r.loc[(plot_data_conc_r['Conc'] > low_cut) & (plot_data_conc_r['Conc'] < high_cut)]
            
            plot_data_l_2 = plot_data_conc_l.loc[plot_data_conc_l['Conc'] >= high_cut]
            plot_data_r_2 = plot_data_conc_r.loc[plot_data_conc_r['Conc'] >= high_cut]
            for i, ax in enumerate(axs.flat):
                
                plot_y = plot_data_conc_l['Integral_trapz']*1e6
                plot_x = plot_data_conc_l['dI_dz']*1e6
                plot_z = plot_data_conc_l['Conc']
                
                vmin=np.min(plot_z)
                vmax=np.max(plot_z)
                
                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ylim_l = ax.get_ylim()
            xlim_l = ax.get_xlim()
            ax.tick_params(direction="in", length=4, width=3)
            fig.colorbar(plot_c, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_001.svg', format='svg', dpi=400, bbox_inches='tight')
            
            Title = 'Left RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            print('Diffusion low conc = {0:.2e}'.format(np.median(plot_data_l_0['D'])))
            for i, ax in enumerate(axs.flat):       
                plot_y_0 = plot_data_l_0['Integral_trapz']*1e6
                plot_x_0 = plot_data_l_0['dI_dz']*1e6
                plot_z_0 = plot_data_l_0['Conc']
                plot_c_0 = ax.scatter(plot_x_0, plot_y_0, c=plot_z_0, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)
                
            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ax.set_ylim(ylim_l[0], ylim_l[1])
            ax.set_xlim(xlim_l[0], xlim_l[1])
            fig.colorbar(plot_c_0, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_002.svg', format='svg', dpi=400, bbox_inches='tight')
                            
            Title = 'Left RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            print('Diffusion mid conc = {0:.2e}'.format(np.median(plot_data_l_1['D'])))
            for i, ax in enumerate(axs.flat):
                plot_y_1 = plot_data_l_1['Integral_trapz']*1e6
                plot_x_1 = plot_data_l_1['dI_dz']*1e6
                plot_z_1 = plot_data_l_1['Conc']
                plot_c_1 = ax.scatter(plot_x_1, plot_y_1, c=plot_z_1, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)
            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ax.set_ylim(ylim_l[0], ylim_l[1])
            ax.set_xlim(xlim_l[0], xlim_l[1])
            fig.colorbar(plot_c_1, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_003.svg', format='svg', dpi=400, bbox_inches='tight')

            Title = 'Left RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            print('Diffusion high conc = {0:.2e}'.format(np.median(plot_data_l_2['D'])))
            for i, ax in enumerate(axs.flat):
                plot_y_2 = plot_data_l_2['Integral_trapz']*1e6
                plot_x_2 = plot_data_l_2['dI_dz']*1e6
                plot_z_2 = plot_data_l_2['Conc']
                plot_c_2 = ax.scatter(plot_x_2, plot_y_2, c=plot_z_2, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ax.set_ylim(ylim_l[0], ylim_l[1])
            ax.set_xlim(xlim_l[0], xlim_l[1])
            fig.colorbar(plot_c_2, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_004.svg', format='svg', dpi=400, bbox_inches='tight')

            'Right: RA vs DF: Conc'
            Title = 'Right: RA vs DF :Conc'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = plot_data_conc_r['Integral_trapz']*1e6
                plot_x = plot_data_conc_r['dI_dz']*1e6
                plot_z = plot_data_conc_r['Conc']
                
                vmin=np.min(plot_z)
                vmax=np.max(plot_z)

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z,cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ylim_r = ax.get_ylim()
            xlim_r = ax.get_xlim()
            fig.colorbar(plot_c, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_005.svg', format='svg', dpi=400, bbox_inches='tight')

            Title = 'Right RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            print('Diffusion low conc = {0:.2e}'.format(np.median(plot_data_r_0['D'])))
            for i, ax in enumerate(axs.flat):       
                plot_y_0 = plot_data_r_0['Integral_trapz']*1e6
                plot_x_0 = plot_data_r_0['dI_dz']*1e6
                plot_z_0 = plot_data_r_0['Conc']
                plot_c_0 = ax.scatter(plot_x_0, plot_y_0, c=plot_z_0, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ax.set_ylim(ylim_r[0], ylim_r[1])
            ax.set_xlim(xlim_r[0], xlim_r[1])
            fig.colorbar(plot_c_0, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_006.svg', format='svg', dpi=400, bbox_inches='tight')

            Title = 'Right RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            print('Diffusion mid conc = {0:.2e}'.format(np.median(plot_data_r_1['D'])))
            for i, ax in enumerate(axs.flat):
                plot_y_1 = plot_data_r_1['Integral_trapz']*1e6
                plot_x_1 = plot_data_r_1['dI_dz']*1e6
                plot_z_1 = plot_data_r_1['Conc']
                plot_c_1 = ax.scatter(plot_x_1, plot_y_1, c=plot_z_1, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ax.set_ylim(ylim_r[0], ylim_r[1])
            ax.set_xlim(xlim_r[0], xlim_r[1])
            fig.colorbar(plot_c_1, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_007.svg', format='svg', dpi=400, bbox_inches='tight')

            Title = 'Right RA vs DF: Conc'
            fig, axs = plt.subplots(squeeze=False)
            print('Diffusion high conc = {0:.2e}'.format(np.median(plot_data_r_2['D'])))
            for i, ax in enumerate(axs.flat):
                plot_y_2 = plot_data_r_2['Integral_trapz']*1e6
                plot_x_2 = plot_data_r_2['dI_dz']*1e6
                plot_z_2 = plot_data_r_2['Conc']
                plot_c_2 = ax.scatter(plot_x_2, plot_y_2, c=plot_z_2, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            ax.set_ylim(ylim_r[0], ylim_r[1])
            ax.set_xlim(xlim_r[0], xlim_r[1])
            fig.colorbar(plot_c_2, ax=ax, label='DNA duplexes/UNC').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_2+str(xtal)+'_008.svg', format='svg', dpi=400, bbox_inches='tight')


        if xtal_fold_i in [4] :
            if xtal == 'C1':
                time_drop = -1
                whisker_outlier_filename_L = 'new_Whisker_outlier/1L_WhiskerPlot_Outliers_Index'
                whisker_outlier_filename_R = 'new_Whisker_outlier/1R_WhiskerPlot_Outliers_Index'
                index_l = pd.read_csv(whisker_outlier_filename_L, delimiter=',')
                index_r = pd.read_csv(whisker_outlier_filename_R, delimiter=',')
                index_l = index_l['0'].tolist()
                index_r = index_r['0'].tolist()   
                index_l = [index for index in index_l if index in range(np.shape(df1l)[0]*(xtal_fold_i-1),np.shape(df1l)[0]*xtal_fold_i)]
                index_r = [index for index in index_r if index in range(np.shape(df1r)[0]*(xtal_fold_i-1),np.shape(df1r)[0]*xtal_fold_i)]
                index_l = [index - np.shape(df1l)[0]*(xtal_fold_i-1)-1 for index in index_l]
                index_r = [index - np.shape(df1r)[0]*(xtal_fold_i-1)-1 for index in index_r]

            elif xtal == 'C2':
                time_drop = -1
                whisker_outlier_filename_L = 'new_Whisker_outlier/2L_WhiskerPlot_Outliers_Index'
                whisker_outlier_filename_R = 'new_Whisker_outlier/2R_WhiskerPlot_Outliers_Index'
                index_l = pd.read_csv(whisker_outlier_filename_L, delimiter=',')
                index_r = pd.read_csv(whisker_outlier_filename_R, delimiter=',')
                index_l = index_l['0'].tolist()
                index_r = index_r['0'].tolist()   
                index_l = [index for index in index_l if index in range(np.shape(df1l)[0]*(xtal_fold_i-1),np.shape(df1l)[0]*xtal_fold_i)]
                index_r = [index for index in index_r if index in range(np.shape(df1r)[0]*(xtal_fold_i-1),np.shape(df1r)[0]*xtal_fold_i)]
                index_l = [index - np.shape(df1l)[0]*(xtal_fold_i-1)-1 for index in index_l]
                index_r = [index - np.shape(df1r)[0]*(xtal_fold_i-1)-1 for index in index_r]
                
            elif xtal == 'C3':
                time_drop = -1
                whisker_outlier_filename_L = 'new_Whisker_outlier/3L_WhiskerPlot_Outliers_Index'
                whisker_outlier_filename_R = 'new_Whisker_outlier/3R_WhiskerPlot_Outliers_Index'
                index_l = pd.read_csv(whisker_outlier_filename_L, delimiter=',')
                index_r = pd.read_csv(whisker_outlier_filename_R, delimiter=',')
                index_l = index_l['0'].tolist()
                index_r = index_r['0'].tolist()   
                index_l = [index for index in index_l if index in range(np.shape(df1l)[0]*(xtal_fold_i-1),np.shape(df1l)[0]*xtal_fold_i)]
                index_r = [index for index in index_r if index in range(np.shape(df1r)[0]*(xtal_fold_i-1),np.shape(df1r)[0]*xtal_fold_i)]
                index_l = [index - np.shape(df1l)[0]*(xtal_fold_i-1)-1 for index in index_l]
                index_r = [index - np.shape(df1r)[0]*(xtal_fold_i-1)-1 for index in index_r]
            
            plot_data_l = tables_data_l
            plot_data_l_03 = tables_data_l.iloc[index_l]
            plot_data_r = tables_data_r
            plot_data_r_03 = tables_data_r.iloc[index_r]
            
            Title = 'Left: D v.s. Conc: Depth'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = plot_data_l['Integral_trapz']*1e6
                plot_x = plot_data_l['dI_dz']*1e6
                plot_z = plot_data_l['depth']
                vmin=np.min(plot_z)
                vmax=np.max(plot_z)

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z,cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)
                
            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            fig.colorbar(plot_c, ax=ax, label='Depth [\u03BCm]').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_1+str(xtal)+'_003.svg', format='svg', dpi=400, bbox_inches='tight')


            Title = 'Right: D v.s. Conc: Depth'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = plot_data_r['Integral_trapz']*1e6
                plot_x = plot_data_r['dI_dz']*1e6
                plot_z = plot_data_r['depth']
                vmin=np.min(plot_z)
                vmax=np.max(plot_z)

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z,cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)

            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            fig.colorbar(plot_c, ax=ax, label='Depth [\u03BCm]').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_1+str(xtal)+'_006.svg', format='svg', dpi=400, bbox_inches='tight')


            'Left + Right: DNA vs z-direction: Depth'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = pd.concat([plot_data_l['Conc'], plot_data_r['Conc'][::-1]])
                plot_x = pd.concat([plot_data_l['depth'],
                    plot_data_r['depth']+np.max(plot_data_r['depth'])])
                plot_z = pd.concat([plot_data_l['depth'],
                    plot_data_r['depth'][::-1]])

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z, cmap=my_cmap, s=10, vmin=vmin, vmax=vmax)

            ax.set_ylabel('DNA duplexes/UCN')
            ax.set_xlabel('z-direction [\u03BCm]')
            ax.tick_params(direction="in", length=4, width=3)
            fig.colorbar(plot_c, ax=ax, label='Depth [\u03BCm]').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_1+str(xtal)+'_001.svg', format='svg', dpi=400, bbox_inches='tight')


            'Left: RA vs DF: Time'
            Title = 'Left: D v.s. Conc: Time'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = plot_data_l['Integral_trapz']*1e6
                plot_x = plot_data_l['dI_dz']*1e6
                plot_z = plot_data_l['time_pt']
                vmin=np.min(plot_z)
                vmax=np.max(plot_z)

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z, cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)
                
            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)
            fig.colorbar(plot_c, ax=ax, label='Time [sec]').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_1+str(xtal)+'_004.svg', format='svg', dpi=400, bbox_inches='tight')


            'Right: RA vs DF: Time'
            Title = 'Right: D v.s. Conc: Time'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = plot_data_r['Integral_trapz']*1e6
                plot_x = plot_data_r['dI_dz']*1e6
                plot_z = plot_data_r['time_pt']
                vmin=np.min(plot_z)
                vmax=np.max(plot_z)

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z,cmap=my_cmap, marker='o', vmin=vmin, vmax=vmax)
                
            ax.set_ylabel('RA [mM·cm/sec]')
            ax.set_xlabel('DF [mM/cm]')
            ax.tick_params(direction="in", length=4, width=3)

            fig.colorbar(plot_c, ax=ax, label='Time [sec]').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_1+str(xtal)+'_007.svg', format='svg', dpi=400, bbox_inches='tight')


            'Left & Right: DNA vs z-direction: Time'
            fig, axs = plt.subplots(squeeze=False)
            for i, ax in enumerate(axs.flat):

                plot_y = pd.concat([plot_data_l['Conc'],
                    plot_data_r['Conc'][::-1]])
                plot_x = pd.concat([plot_data_l['depth'],
                    plot_data_r['depth']+np.max(plot_data_r['depth'])])
                plot_z = pd.concat([plot_data_l['time_pt'],
                    plot_data_r['time_pt'][::-1]])

                length = len(plot_z)
                my_cmap = cmap_generation(length)
                plot_c = ax.scatter(plot_x, plot_y, c=plot_z, cmap=my_cmap, s=10, vmin=vmin, vmax=vmax)
                
            ax.set_ylabel('DNA duplexes/UCN')
            ax.set_xlabel('z-direction [\u03BCm]')
            ax.tick_params(direction="in", length=4, width=3)
            fig.colorbar(plot_c, ax=ax, label='Time [sec]').ax.tick_params(
                length=4, width=3)
            if xtal_fold_i == 4:
                plt.show()
                fig.savefig(save_file_path_1+str(xtal)+'_002.svg', format='svg', dpi=400, bbox_inches='tight')


            if xtal == 'C1':
                Title = 'Anomalous Left'
                fig = plt.figure(figsize=(3, 4))
                ax = fig.add_axes([0, 0, 1, 1])

                plot_data_1_1_l = tables_data_l.loc[tables_data_l['time_pt'] <= t_eval[2]]
                ax.plot(plot_data_l['dI_dz']*1e6,
                        plot_data_l['Integral_trapz']*1e6, 'o', color='darkgray',  fillstyle='none')
                plot_data_l_1 = tables_data_l.loc[tables_data_l['time_pt'] == t_eval[-2]]
                ax.plot(plot_data_l_1['dI_dz']*1e6,
                        plot_data_l_1['Integral_trapz']*1e6, 'bo', fillstyle='none')
                
                plot_data_l_1 = tables_data_l.loc[tables_data_l['time_pt'] == t_eval[-3]]
                ax.plot(plot_data_l_1['dI_dz']*1e6,
                        plot_data_l_1['Integral_trapz']*1e6, 'go', fillstyle='none')
                
                plot_data_l_1 = tables_data_l.loc[tables_data_l['time_pt'] == t_eval[-4]]
                ax.plot(plot_data_l_1['dI_dz']*1e6,
                        plot_data_l_1['Integral_trapz']*1e6, 'yo', fillstyle='none')
                if save_fig_3 == True:
                    ax.plot(plot_data_l_03['dI_dz']*1e6,
                                plot_data_l_03['Integral_trapz']*1e6, 'mo', label='Whisker Outliers')
                if save_fig_2 == True:
                    ax.plot([-5, 20], [0.0e-8, 8e-8], 'r--', label='Max {0:.2e}'.format((8e-8-0e-8)/(20-(-5))))
                    ax.plot([2, 37], [0.0e-8, 6e-8], 'r--', label='Min {0:.2e}'.format((6e-8-0e-8)/35))
                plt.tick_params(left=False, right=False, labelleft=False,
                                labelbottom=False, bottom=False)
                plt.title(Title)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(
                    1, 1.05), shadow=True, mode=None, fancybox=True)
                if xtal_fold_i == 4:
                    plt.show()
                    fig.savefig(save_file_path_1+str(xtal)+'_005.svg', format='svg', dpi=400, bbox_inches='tight')


                Title = 'Anomalous Right'
                fig = plt.figure(figsize=(3, 4))
                ax = fig.add_axes([0, 0, 1, 1])
                plot_data_1_r = tables_data_r.drop(tables_data_r.loc[tables_data_r['time_pt'] == t_eval[-4]].index)

                plot_data_1_2_r = tables_data_r.loc[tables_data_r['time_pt'] == t_eval[-4]]
                ax.plot(plot_data_r['dI_dz']*1e6,
                        plot_data_r['Integral_trapz']*1e6, 'o', color='darkgray',  fillstyle='none')
                plot_data_r_1 = tables_data_r.loc[tables_data_r['time_pt'] == t_eval[-2]]
                ax.plot(plot_data_r_1['dI_dz']*1e6,
                        plot_data_r_1['Integral_trapz']*1e6, 'bo', fillstyle='none')
                
                plot_data_r_1 = tables_data_r.loc[tables_data_r['time_pt'] == t_eval[-3]]
                ax.plot(plot_data_r_1['dI_dz']*1e6,
                        plot_data_r_1['Integral_trapz']*1e6, 'go', fillstyle='none')
                
                plot_data_r_1 = tables_data_r.loc[tables_data_r['time_pt'] == t_eval[-4]]
                ax.plot(plot_data_r_1['dI_dz']*1e6,
                        plot_data_r_1['Integral_trapz']*1e6, 'yo', fillstyle='none')
                if save_fig_3 == True:
                    ax.plot(plot_data_r_03['dI_dz']*1e6,
                                plot_data_r_03['Integral_trapz']*1e6, 'mo', label='Whisker Outliers')
                if save_fig_2 == True:
                    ax.plot([-2, 28], [0e-8, 8e-8], 'r--', label='Max {0:.2e}'.format((8e-8-0e-8)/(28-(-2))))
                    ax.plot([9, 49], [0e-8, 6e-8], 'r--', label='Min {0:.2e}'.format((6e-8-0e-8)/(49-9)))
                plt.tick_params(left=False, right=False, labelleft=False,
                                labelbottom=False, bottom=False)
                plt.title(Title)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(
                    1, 1.05), shadow=True, mode=None, fancybox=True)
                if xtal_fold_i == 4:
                    plt.show()
                    fig.savefig(save_file_path_1+str(xtal)+'_008.svg', format='svg', dpi=400, bbox_inches='tight')


            elif xtal == 'C2':
                Title = 'Anomalous Left'
                fig = plt.figure(figsize=(3, 4))
                ax = fig.add_axes([0, 0, 1, 1])
                ax.plot(plot_data_l['dI_dz']*1e6,
                        plot_data_l['Integral_trapz']*1e6, 'o', color='darkgray', fillstyle='none')
                if save_fig_3 == True:
                    ax.plot(plot_data_l_03['dI_dz']*1e6,
                                plot_data_l_03['Integral_trapz']*1e6, 'mo', label='Whisker Outliers')
                if save_fig_2 == True:
                    ax.plot([-1, 4], [0.0e-8, 2.5e-8], 'r--', label='Max {0:.2e}'.format((2.5e-8-0e-8)/5))
                    ax.plot([1, 17], [0.0e-8, 2.5e-8], 'r--', label='Min {0:.2e}'.format((2.5e-8-0.0e-8)/16))
                plt.tick_params(left=False, right=False, labelleft=False,
                                labelbottom=False, bottom=False)
                plt.title(Title)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(
                    1, 1.05), shadow=True, mode=None, fancybox=True)
                if xtal_fold_i == 4:
                    plt.show()
                    fig.savefig(save_file_path_1+str(xtal)+'_005.svg', format='svg', dpi=400, bbox_inches='tight')


                Title = 'Anomalous Right'
                fig = plt.figure(figsize=(3, 4))
                ax = fig.add_axes([0, 0, 1, 1])
                plot_data_2 = tables_data_r
                ax.plot(plot_data_r['dI_dz']*1e6,
                        plot_data_r['Integral_trapz']*1e6, 'o', color='darkgray', fillstyle='none')
                if save_fig_3 == True:
                    ax.plot(plot_data_r_03['dI_dz']*1e6,
                                plot_data_r_03['Integral_trapz']*1e6, 'mo', label='Whisker Outliers')
                if save_fig_2 == True:
                    ax.plot([-1, 5], [0.0e-8, 2.5e-8], 'r--', label='Max {0:.2e}'.format((2.5e-8-0e-8)/6))
                    ax.plot([1, 18], [0.0e-8, 2.5e-8], 'r--', label='Min {0:.2e}'.format((2.5e-8-0.0e-8)/17))
                plt.tick_params(left=False, right=False, labelleft=False,
                                labelbottom=False, bottom=False)
                plt.title(Title)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(
                    1, 1.05), shadow=True, mode=None, fancybox=True)
                if xtal_fold_i == 4:
                    plt.show()
                    fig.savefig(save_file_path_1+str(xtal)+'_008.svg', format='svg', dpi=400, bbox_inches='tight')


            elif xtal == 'C3':
                Title = 'Anomalous Left'
                fig = plt.figure(figsize=(3, 4))
                ax = fig.add_axes([0, 0, 1, 1])
                ax.plot(plot_data_l['dI_dz']*1e6,
                        plot_data_l['Integral_trapz']*1e6, 'o', color='darkgray', fillstyle='none')
                plot_data_l_1 = tables_data_l.loc[tables_data_l['time_pt'] == t_eval[-2]]
                ax.plot(plot_data_l_1['dI_dz']*1e6,
                        plot_data_l_1['Integral_trapz']*1e6, 'bo', fillstyle='none')
                if save_fig_3 == True:
                    ax.plot(plot_data_l_03['dI_dz']*1e6,
                            plot_data_l_03['Integral_trapz']*1e6, 'mo', label='Whisker Outliers')
                if save_fig_2 == True:
                    ax.plot([0, 11], [0.0e-7, 2e-7], 'r--', label='Max {0:.2e}'.format((2e-7-0e-7)/11))
                    ax.plot([0, 23], [0.0e-7, 1.1e-7], 'r--', label='Min {0:.2e}'.format((1.1e-7-0e-7)/23))
                plt.tick_params(left=False, right=False, labelleft=False,
                                labelbottom=False, bottom=False)
                plt.title(Title)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(
                    1, 1.05), shadow=True, mode=None, fancybox=True)
                if xtal_fold_i == 4:
                    plt.show()
                    fig.savefig(save_file_path_1+str(xtal)+'_005.svg', format='svg', dpi=400, bbox_inches='tight')


                Title = 'Anomalous Right'
                fig = plt.figure(figsize=(3, 4))
                ax = fig.add_axes([0, 0, 1, 1])
                plot_data_3 = tables_data_r
                ax.plot(plot_data_r['dI_dz']*1e6,
                        plot_data_r['Integral_trapz']*1e6, 'o', color='darkgray', fillstyle='none')
                plot_data_r_1 = tables_data_r.loc[tables_data_r['time_pt'] == t_eval[-2]]
                ax.plot(plot_data_r_1['dI_dz']*1e6,
                        plot_data_r_1['Integral_trapz']*1e6, 'bo', fillstyle='none')
                if save_fig_3 == True:
                    ax.plot(plot_data_l_03['dI_dz']*1e6,
                                plot_data_l_03['Integral_trapz']*1e6, 'mo', label='Whisker Outliers')
                if save_fig_2 == True:
                    ax.plot([-0.5, 10.5], [0.0e-7, 2.2e-7], 'r--', label='Max {0:.2e}'.format((2.2e-7-0e-7)/11))
                    ax.plot([1.5, 19.5], [0.0e-7, 1.4e-7], 'r--', label='Min {0:.2e}'.format((1.4e-7-0e-7)/18))
                plt.tick_params(left=False, right=False, labelleft=False,
                                labelbottom=False, bottom=False)
                plt.title(Title)
                lgd = plt.legend(loc='upper left', bbox_to_anchor=(
                    1, 1.05), shadow=True, mode=None, fancybox=True)
                if xtal_fold_i == 4:
                    plt.show()
                    fig.savefig(save_file_path_1+str(xtal)+'_008.svg', format='svg', dpi=400, bbox_inches='tight')

    df1l.to_csv(filename1_l, index=False)
    df2l.to_csv(filename2_l, index=False)
    df3l.to_csv(filename3_l, index=False)
    df4l.to_csv(filename4_l, index=False)

    df1r.to_csv(filename1_r, index=False)
    df2r.to_csv(filename2_r, index=False)
    df3r.to_csv(filename3_r, index=False)
    df4r.to_csv(filename4_r, index=False)
print('Finished.')
# %%
