# %%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl

'---Image settings-----'
size = 25
plt.rc('axes', titlesize=size)     # fontsize of the axes title
plt.rc('axes', labelsize=size)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=size)    # fontsize of the tick labels
plt.rc('ytick', labelsize=size)    # fontsize of the tick labels
plt.rc('legend', fontsize=size)    # legend fontsize
plt.rc('figure', titlesize=5)      # fontsize of the figure title
mpl.rcParams['figure.figsize'] = 4, 5
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10

file_path = 'Whisker_plot_csv/'

filename1_l = file_path+'all_crystal_20fold_smooth_c1_Dl.csv'
filename1t_l = file_path+'all_crystal_20fold_smooth_c1_Timel.csv'
filename1c_l = file_path+'all_crystal_20fold_smooth_c1_Concl.csv'
filename1d_l = file_path+'all_crystal_20fold_smooth_c1_Depthl.csv'

filename1_r = file_path+'all_crystal_20fold_smooth_c1_Dr.csv'
filename1t_r = file_path+'all_crystal_20fold_smooth_c1_Timer.csv'
filename1c_r = file_path+'all_crystal_20fold_smooth_c1_Concr.csv'
filename1d_r = file_path+'all_crystal_20fold_smooth_c1_Depthr.csv'

data1_l = pd.read_csv(filename1_l, delimiter=',')
data1t_l = pd.read_csv(filename1t_l, delimiter=',')
data1c_l = pd.read_csv(filename1c_l, delimiter=',')
data1d_l = pd.read_csv(filename1d_l, delimiter=',')

data1_r = pd.read_csv(filename1_r, delimiter=',')
data1t_r = pd.read_csv(filename1t_r, delimiter=',')
data1c_r = pd.read_csv(filename1c_r, delimiter=',')
data1d_r = pd.read_csv(filename1d_r, delimiter=',')

filename2_l = file_path+'all_crystal_20fold_smooth_c2_Dl.csv'
filename2t_l = file_path+'all_crystal_20fold_smooth_c2_Timel.csv'
filename2c_l = file_path+'all_crystal_20fold_smooth_c2_Concl.csv'
filename2d_l = file_path+'all_crystal_20fold_smooth_c2_Depthl.csv'

filename2_r = file_path+'all_crystal_20fold_smooth_c2_Dr.csv'
filename2t_r = file_path+'all_crystal_20fold_smooth_c2_Timer.csv'
filename2c_r = file_path+'all_crystal_20fold_smooth_c2_Concr.csv'
filename2d_r = file_path+'all_crystal_20fold_smooth_c2_Depthr.csv'

data2_l = pd.read_csv(filename2_l, delimiter=',')
data2t_l = pd.read_csv(filename2t_l, delimiter=',')
data2c_l = pd.read_csv(filename2c_l, delimiter=',')
data2d_l = pd.read_csv(filename2d_l, delimiter=',')

data2_r = pd.read_csv(filename2_r, delimiter=',')
data2t_r = pd.read_csv(filename2t_r, delimiter=',')
data2c_r = pd.read_csv(filename2c_r, delimiter=',')
data2d_r = pd.read_csv(filename2d_r, delimiter=',')

filename3_l = file_path+'all_crystal_20fold_smooth_c3_Dl.csv'
filename3t_l = file_path+'all_crystal_20fold_smooth_c3_Timel.csv'
filename3c_l = file_path+'all_crystal_20fold_smooth_c3_Concl.csv'
filename3d_l = file_path+'all_crystal_20fold_smooth_c3_Depthl.csv'

filename3_r = file_path+'all_crystal_20fold_smooth_c3_Dr.csv'
filename3t_r = file_path+'all_crystal_20fold_smooth_c3_Timer.csv'
filename3c_r = file_path+'all_crystal_20fold_smooth_c3_Concr.csv'
filename3d_r = file_path+'all_crystal_20fold_smooth_c3_Depthr.csv'

data3_l = pd.read_csv(filename3_l, delimiter=',')
data3t_l = pd.read_csv(filename3t_l, delimiter=',')
data3c_l = pd.read_csv(filename3c_l, delimiter=',')
data3d_l = pd.read_csv(filename3d_l, delimiter=',')

data3_r = pd.read_csv(filename3_r, delimiter=',')
data3t_r = pd.read_csv(filename3t_r, delimiter=',')
data3c_r = pd.read_csv(filename3c_r, delimiter=',')
data3d_r = pd.read_csv(filename3d_r, delimiter=',')

'Crystal 1 Left'
data1__l = data1_l.iloc[:, 0]
data1t__l = data1t_l.iloc[:, 0]
data1c__l = data1c_l.iloc[:, 0]
data1d__l = data1d_l.iloc[:, 0]
med_data1_l = [np.median(data1__l)]
for j in range(7):
    k = j+1
    data1__l = pd.concat([data1__l, data1_l.iloc[:, k]], ignore_index=True)
    data1t__l = pd.concat([data1t__l, data1t_l.iloc[:, k]], ignore_index=True)
    data1c__l = pd.concat([data1c__l, data1c_l.iloc[:, k]], ignore_index=True)
    data1d__l = pd.concat([data1d__l, data1d_l.iloc[:, k]], ignore_index=True)
    med_data1_l.append(np.median(data1_l.iloc[:, k]))
print('Crystal 1 Left')
print(med_data1_l)
print('Median: {:.2e}'.format(np.median(data1__l)))

'Crystal 1 Right'
data1__r = data1_r.iloc[:, 0]
data1t__r = data1t_r.iloc[:, 0]
data1c__r = data1c_r.iloc[:, 0]
data1d__r = data1d_r.iloc[:, 0]
med_data1_r = [np.median(data1__r)]
for j in range(7):
    k = j+1
    data1__r = pd.concat([data1__r, data1_r.iloc[:, k]], ignore_index=True)
    data1t__r = pd.concat([data1t__r, data1t_r.iloc[:, k]], ignore_index=True)
    data1c__r = pd.concat([data1c__r, data1c_r.iloc[:, k]], ignore_index=True)
    data1d__r = pd.concat([data1d__r, data1d_r.iloc[:, k]], ignore_index=True)
    med_data1_r.append(np.median(data1_r.iloc[:, k]))
print('Crystal 1 Right')
print(med_data1_r)
print('Median: {:.2e}'.format(np.median(data1__r)))

mean = np.mean([med_data1_l, med_data1_r])
std = np.std([med_data1_l, med_data1_r])
print('Crystal 1 diffusion rate = {0:.1e}±{1:.1e}\n'.format(mean, std))

'Crystal 2 Left'
data2__l = data2_l.iloc[:, 0]
data2t__l = data2t_l.iloc[:, 0]
data2c__l = data2c_l.iloc[:, 0]
data2d__l = data2d_l.iloc[:, 0]
med_data2_l = [np.median(data2__l)]
for j in range(7):
    k = j+1
    data2__l = pd.concat([data2__l, data2_l.iloc[:, k]], ignore_index=True)
    data2t__l = pd.concat([data2t__l, data2t_l.iloc[:, k]], ignore_index=True)
    data2c__l = pd.concat([data2c__l, data2c_l.iloc[:, k]], ignore_index=True)
    data2d__l = pd.concat([data2d__l, data2d_l.iloc[:, k]], ignore_index=True)
    med_data2_l.append(np.median(data2_l.iloc[:, k]))
print('Crystal 2 Left')
print(med_data2_l)
print('Median: {:.2e}'.format(np.median(data2__l)))

'Crystal 2 Right'
data2__r = data2_r.iloc[:, 0]
data2t__r = data2t_r.iloc[:, 0]
data2c__r = data2c_r.iloc[:, 0]
data2d__r = data2d_r.iloc[:, 0]
med_data2_r = [np.median(data2__r)]
for j in range(7):
    k = j+1
    data2__r = pd.concat([data2__r, data2_r.iloc[:, k]], ignore_index=True)
    data2t__r = pd.concat([data2t__r, data2t_r.iloc[:, k]], ignore_index=True)
    data2c__r = pd.concat([data2c__r, data2c_r.iloc[:, k]], ignore_index=True)
    data2d__r = pd.concat([data2d__r, data2d_r.iloc[:, k]], ignore_index=True)
    med_data2_r.append(np.median(data2_r.iloc[:, k]))
print('Crystal 2 Right')
print(med_data2_r)
print('Median: {:.2e}'.format(np.median(data2__r)))

mean = np.mean([med_data2_l, med_data2_r])
std = np.std([med_data2_l, med_data2_r])
print('Crystal 2 diffusion rate = {0:.1e}±{1:.1e}\n'.format(mean, std))

'Crystal 3 Left'
data3__l = data3_l.iloc[:, 0]
data3t__l = data3t_l.iloc[:, 0]
data3c__l = data3c_l.iloc[:, 0]
data3d__l = data3d_l.iloc[:, 0]
med_data3_l = [np.median(data3__l)]
for j in range(7):
    k = j+1
    data3__l = pd.concat([data3__l, data3_l.iloc[:, k]], ignore_index=True)
    data3t__l = pd.concat([data3t__l, data3t_l.iloc[:, k]], ignore_index=True)
    data3c__l = pd.concat([data3c__l, data3c_l.iloc[:, k]], ignore_index=True)
    data3d__l = pd.concat([data3d__l, data3d_l.iloc[:, k]], ignore_index=True)
    med_data3_l.append(np.median(data3_l.iloc[:, k]))
print('Crystal 3 Left')
print(med_data3_l)
print('Median: {:.2e}'.format(np.median(data3__l)))

'Crystal 3 Right'
data3__r = data3_r.iloc[:, 0]
data3t__r = data3t_r.iloc[:, 0]
data3c__r = data3c_r.iloc[:, 0]
data3d__r = data3d_r.iloc[:, 0]
med_data3_r = [np.median(data3__r)]
for j in range(7):
    k = j+1
    data3__r = pd.concat([data3__r, data3_r.iloc[:, k]], ignore_index=True)
    data3t__r = pd.concat([data3t__r, data3t_r.iloc[:, k]], ignore_index=True)
    data3c__r = pd.concat([data3c__r, data3c_r.iloc[:, k]], ignore_index=True)
    data3d__r = pd.concat([data3d__r, data3d_r.iloc[:, k]], ignore_index=True)
    med_data3_r.append(np.median(data3_r.iloc[:, k]))
print('Crystal 3 Right')
print(med_data3_r)
print('Median: {:.2e}'.format(np.median(data3__r)))
mean = np.mean([med_data3_l, med_data3_r])
std = np.std([med_data3_l, med_data3_r])
print('Crystal 3 diffusion rate = {0:.1e}±{1:.1e}\n'.format(mean, std))

fig = plt.figure()
a1 = fig.add_axes([0, 0, 1, 1])
plot_data = [data1__l, data1__r, data2__l, data2__r, data3__l, data3__r]
box = a1.boxplot(plot_data, showfliers=True)
plt.ylabel('$D_0$'+' ['+'$cm^2$'+'/sec]')
my_xticks = ['1L', '1R', '2L', '2R', '3L', '3R']
plt.xticks([1, 2, 3, 4, 5, 6], my_xticks)
plt.xlabel('Crystal #')
plt.tick_params(direction="in", length=4, width=3)
plt.show()

plt.ion()
mpl.is_interactive()
fig = plt.figure()
a1 = fig.add_axes([0, 0, 1, 1])
plot_data = [data1__l, data1__r, data2__l, data2__r, data3__l, data3__r]
a1.boxplot(plot_data, showfliers=False)
plt.ylabel('$D_0$'+' ['+'$cm^2$'+'/sec]')
my_xticks = ['1L', '1R', '2L', '2R', '3L', '3R']
plt.xticks([1, 2, 3, 4, 5, 6], my_xticks)
plt.xlabel('Crystal #')
plt.tick_params(direction="in", length=4, width=3)
plt.tight_layout()
plt.show()
# %%
