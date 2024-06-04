# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 07:55:52 2021
@author: Ashlyn
"""
import numpy as np
from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import make_interp_spline, UnivariateSpline
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image
from skimage import io
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter, zoom
import csv

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
# 
' ======================================================================= '
' --- User Defined Settings --------------------------------------------- '
' ==========================/============================================= '

' - - STRING VARIABLES  - - - - - - - - - - - - - - - - - - - - - - - - - '
file_path = 'Crystals_images/'
save_file_path = 'savefig/01_Processed_Image/'

xtal = 'C1'     
file_name = 'Crystal1'

# xtal = 'C2'
# file_name = 'Crystal2'

# xtal = 'C3'
# file_name = 'Crystal3'


save_path_1 = 'new_Crystal_smooth_csv/' + file_name
save_path_1_l = 'new_Crystal_smooth_csv/' + file_name  # File path to left boundary images csv data
save_path_1_r = 'new_Crystal_smooth_csv/' + file_name  # File path to right boundary images csv data
save_path_2 = 'new_Crystal_csv/' + file_name           # File path to images csv data
save_path_2_l = 'new_Crystal_csv/' + file_name # File path to left boundary images csv data
save_path_2_r = 'new_Crystal_csv/' + file_name # File path to right boundary images csv data

fnameTiff = file_path + file_name + '.tiff'   # Multidimensional Tiff data file
fnameRGB = file_path + file_name + '.tif'     # RBG Tiff file for display
' - - NUMERIC VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - '
el_num = 100                # Constrain fit to this number of elements
el_time = 16
pixel_length = 0.68         # um/pixel for JDS's crystal data
conc_conversion = 845351.18 # mole/mL --> DNA duplexes/UCN (unit cell nanopore)

'Setting for each crystal'
if xtal == 'C1':
    trimmed_pixel = 30   # Pixel of image to trim out boundary optical effect
    height = 80          # Height of data slice (pixel)
    t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
              1087, 1207, 1327, 1447, 1567, 1687, 1807]   # Evaluation time points
    left = 170           # Start point of left side of the crystal (pixel)
    right = 270          # Start point of right side of the crystal (pixel)
    center = 186         # Center of the crystal (pixel)
    z_position1 = [1, 5, 19]
    z_position2 = [4, 14, 49]
    xlim0 = [-2, 16]
    xlim1 = [-2, 67]     # x-axis range for fig0 (pixel)
    xlim2 = [-1, 47]     # x-axis range for fig1 (pixel)
    ylim1 = [7000, 55000]
    ylim2 = [0, 8e-8]
    image0 = file_path+'Crystal1_0000.jpg'  # image path for plotting
    image1 = file_path+'Crystal1_0008.jpg'  # image path for plotting
    image2 = file_path+'Crystal1_0015.jpg'  # image path for plotting
    t_num = [0, 8, 15]   # Time index for image0, image1, and image 2
    trim_right = 152     # Pixel of right side of image to trim out for plotting
    trim_left = 152      # Pixel of left side of image to trim out for plotting
    trim_top = 100       # Pixel of top side of image to trim out for plotting
    trim_bottom = 120    # Pixel of bottom side of image to trim out for plotting
    sigma = (1, 2)

elif xtal == 'C2':
    trimmed_pixel = 40  # Pixel of image to trim out boundary optical effect
    height = 100        # Height of data slice (pixel)
    t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
              1087, 1207, 1327, 1447, 1567, 1687, 1807]  # Evaluation time points
    left = 140          # Start point of left side of the crystal (pixel)
    right = 250         # Start point of right side of the crystal (pixel)
    center = 210        # Center of the crystal (pixel)
    z_position1 = [1, 4, 14]
    z_position2 = [4, 14, 49]
    xlim0 = [-2, 16]
    xlim1 = [-2, 53]    # x-axis range for fig0 (pixel)
    xlim2 = [-1, 36]    # x-axis range for fig1 (pixel)
    ylim1 = [5000, 40000]
    ylim2 = [0, 5e-8]
    image0 = file_path+'Crystal2_0000.jpg'  # image path for plotting
    image1 = file_path+'Crystal2_0008.jpg'  # image path for plotting
    image2 = file_path+'Crystal2_0015.jpg'  # image path for plotting
    t_num = [0, 8, 15]  # Time index for image0, image1, and image 2
    trim_right = 123    # Pixel of right side of image to trim out for plotting
    trim_left = 123     # Pixel of left side of image to trim out for plotting
    trim_top = 90       # Pixel of top side of image to trim out for plotting
    trim_bottom = 70    # Pixel of bottom side of image to trim out for plotting
    sigma = (1, 1)

elif xtal == 'C3':
    trimmed_pixel = 15  # Pixel of image to trim out boundary optical effect
    height = 45         # Height of data slice             
    t_eval = [17, 137, 257, 377, 497, 617, 737, 857, 977,
              1097, 1217, 1337, 1457, 1577, 1697, 1817]    # Evaluation time points
    left = 83           # Start point of left side of the crystal (pixel)
    right = 204         # Start point of right side of the crystal (pixel)
    center = 190        # Center of the crystal (pixel)
    z_position1 = [4, 13, 45]
    z_position2 = [4, 14, 49]
    xlim0 = [-2, 16]
    xlim1 = [-3, 155]   # x-axis range for fig0 (pixel)
    xlim2 = [-1, 110]   # x-axis range for fig1 (pixel)
    ylim1 = [7000, 55000]
    ylim2 = [0, 1e-7]
    image0 = file_path+'Crystal3_0000.jpg' # image path for plotting
    image1 = file_path+'Crystal3_0008.jpg' # image path for plotting
    image2 = file_path+'Crystal3_0015.jpg' # image path for plotting
    t_num = [0, 8, 15]  # Time index for image0, image1, and image 2
    trim_right = 70     # Pixel of right side of image to trim out for plotting
    trim_left = 70      # Pixel of left side of image to trim out for plotting
    trim_top = 60       # Pixel of top side of image to trim out for plotting
    trim_bottom = 50    # Pixel of bottom side of image to trim out for plotting
    sigma = (1, 0.1)

' ======================================================================= '
' --- Gather Image Data ------------------------------------------------- '
' ======================================================================= '

' Determine image dimensions'
columns, rows = Image.open(fnameRGB).size
Img_tiff = Image.open(fnameTiff)
num_images = Img_tiff.n_frames
print('number of images = {0}'.format(num_images))

'These are a guess :(...'
L0 = np.abs(right - left)
print('Crystal original length (um) = {0:.2f}'.format(L0*pixel_length))

xtal_l = left + trimmed_pixel
xtal_r = right - trimmed_pixel
L1 = xtal_r - xtal_l
print('Crystal trimmed length (um) = {0:.2f}'.format(L1*pixel_length))

top = round(center+height)
bottom = round(center-height)

' Gather all image data based on edge deffinitions '
print('Gathering image data...')
im_array = np.zeros([num_images, top-bottom, columns])
xtal_I_vec = np.zeros([num_images, xtal_r-xtal_l])
full_I_vec = np.zeros([num_images, columns])
image = io.imread(fnameTiff)

i = -1
for image_path in [image0, image1, image2]:
    i += 1
    jpg_image = io.imread(image_path)

    fig = plt.figure()
    a1 = fig.add_axes([0, 0, 1, 1])
    a1.plot([xtal_l-trim_left, xtal_l-trim_left], [center+height-trim_top,
            center-height-trim_top], 'w--', linewidth=3)
    a1.plot([xtal_r-trim_right, xtal_r-trim_right], [center+height-trim_top,
            center-height-trim_top], 'w--', linewidth=3)
    a1.plot([(xtal_l-trim_left+xtal_r-trim_right)/2, (xtal_l-trim_left+xtal_r-trim_right)/2], [center+height-trim_top,
            center-height-trim_top], 'w--', linewidth=3)
    for fold in np.arange(11):
        a1.plot([xtal_l-trim_left, xtal_r-trim_right], [((bottom-trim_top-1)+int(fold*2*height/10)), ((bottom-trim_top-1)+int(fold*2*height/10))],
                'w--', linewidth=3)
    plt.imshow(jpg_image[trim_top:-trim_bottom, trim_left:-trim_right])

    plt.axis('off')
    position1, position2, channel = np.shape(
        jpg_image[trim_top:-trim_bottom, trim_left:-trim_right])
    a1.text(position2*0.98, position1*0.98,  '%.0f sec' %
            (t_eval[t_num[i]]), size=25, color='white', ha='right')
    plt.show()

image = io.imread(fnameTiff)
I0 = np.min(image[0, bottom-1:top-1, xtal_l:xtal_r])

for fold in np.arange(10):
    for inum in range(0, num_images):
        cur_image = image[inum, :, :]
        im_array[inum, :, :] = cur_image[bottom-1:top-1, :]
        full_I_vec[inum, :] = np.mean(
            cur_image[((bottom-1)+int(fold*2*height/10)):((bottom-1)+int((fold+1)*2*height/10)), :], axis=0)
        xtal_I_vec[inum, :] = full_I_vec[inum, xtal_l:xtal_r]
    
    print('xtal_I_vec=',np.shape(xtal_I_vec))
    if fold == 4:
        t_axis = np.linspace(0, t_eval[-1], el_time, endpoint=True)  # cm -> um
        z_axis = np.linspace(0, Lp, np.shape(xtal_I_vec)[1], endpoint=True)  # cm -> um


        Title = file_name
        fig = plt.figure()
        a0 = fig.add_axes([0, 0, 1, 1])
        lines = []
        annotations = []
        for zi in z_position1: #0-39
            lines.append(a0.plot(xtal_I_vec[:, zi], '-', label="{0:.1f} μm".format(z_axis[i]*10**4), linewidth=3, markersize=12)[0])
            annotations.append("{0:.1f} μm".format(z_axis[zi]*10**4))

        idx = 12
        for line, label in zip(lines, annotations):
            idx -= 1
            # Get y-value of the line at the closest x-value
            y_pos = line.get_ydata()[idx]
            a0.annotate(label, xy=(xlim0[0]+4, y_pos), xytext=(xlim0[0]+4, y_pos), color=line.get_color(),
                        fontsize=25, ha='center')

        plt.xlabel('pixel #')
        plt.ylabel('Intensity')
        plt.ylim([ylim1[0], ylim1[1]])
        plt.xlim([xlim0[0], xlim0[1]])        
        plt.tick_params(direction="in", length=4, width=2)
        plt.show()
        fig.savefig(save_file_path+str(xtal)+'_001.svg', format='svg', dpi=400, bbox_inches='tight')

        fig = plt.figure()
        a0 = fig.add_axes([0, 0, 1, 1])
        lines = []
        annotations = []
        for ti in [0, 3, 7, 11, 15]:
            lines.append(a0.plot(xtal_I_vec[ti, :], '-', label=f"t={t_eval[ti]} sec", linewidth=3, markersize=12)[0])
            annotations.append(f"{t_eval[ti]} sec")

        x_pos = xlim1[1]
        for line, label in zip(lines, annotations):
            x = np.arange(0, np.shape(xtal_I_vec[ti,:])[0], 1)
            # Index of closest x-value to x_pos
            idx = np.abs(x - x_pos).argmin()
            # Get y-value of the line at the closest x-value
            y_pos = line.get_ydata()[idx]
            a0.annotate(label, xy=(xlim1[1]-(xlim1[1]-x[idx])/2, y_pos), xytext=(xlim1[1]-(xlim1[1]-x[idx])/2, y_pos), color=line.get_color(),
                        fontsize=25, ha='center')

        plt.xlabel('pixel #')
        plt.ylabel('Intensity')
        plt.ylim([ylim1[0], ylim1[1]])
        plt.xlim([xlim1[0], xlim1[1]])
        plt.tick_params(direction="in", length=4, width=2)
        plt.show()
        fig.savefig(save_file_path+str(xtal)+'_003.svg', format='svg', dpi=400, bbox_inches='tight')

    'Convert from intensity values to concentrations based on standard'
    xtal_I = xtal_I_vec-I0
    xtal_I_vec_conc = 3.63042E-17*xtal_I**2 + 1.24742E-13*xtal_I + 3.60000E-10

    ' ======================================================================= '
    ' --- Fit Image Data to Diffusion Coef. --------------------------------- '
    ' ======================================================================= '
    Lp = round((xtal_r-xtal_l))*1e-4*pixel_length
    
    
    
    # Apply Gaussian filter
    smoothed_data = gaussian_filter(xtal_I_vec_conc, sigma=sigma)
    # Calculate the zoom factor for each dimension
    zoom_factor = (el_time/np.shape(smoothed_data)[0], 1)
    expanded_data = zoom(smoothed_data, zoom_factor)
    
    zoom_factor = (el_time/np.shape(smoothed_data)[0], el_num/np.shape(smoothed_data)[1])
    spline_dif_mat = zoom(xtal_I_vec_conc, zoom_factor)
    
    dif_mat = expanded_data
    x_dist = np.shape(dif_mat)[1]

    ' Matrix of concentrations in time vs position. '
    new_dif_mat = np.zeros([el_time, el_num])
    new_dif_mat_smooth = np.zeros([el_time, el_num])
    V0 = np.arange(1, x_dist+1) 
    V1 = np.linspace(1, x_dist, el_num, endpoint=True)

    for inum in range(0, el_time):
        c_i_v = dif_mat[inum, :]
        ' Interpolate to find values of c_i_v at coordinates V1.'
        # B-spline degree. Default is cubic, k=3.
        Spline = make_interp_spline(
            V0, uniform_filter1d(c_i_v, size=1), k=3)

        new_dif_mat[inum, :] = Spline(V1)

        Spline = UnivariateSpline(
            V0, uniform_filter1d(c_i_v, size=1), k=3)
        new_dif_mat_smooth[inum, :] = Spline(V1)

    smooth_dif_mat = new_dif_mat_smooth

    'Smooth for each time'
    dif_mat = smooth_dif_mat
    x_dist = np.shape(dif_mat)[1]

    ' Matrix of concentrations in time vs position. '
    new_dif_mat_2 = np.zeros([el_time, el_num])
    new_dif_mat_smooth_2 = np.zeros([el_time, el_num])
    V0 = np.arange(1, x_dist+1) 
    V1 = np.linspace(1, x_dist, el_num, endpoint=True)

    for inum in range(0, el_time):
        c_i_v = dif_mat[inum, :]
        ' Interpolate to find values of c_i_v at coordinates V1.'
        # B-spline degree. Default is cubic, k=3.
        Spline = make_interp_spline(
            V0, uniform_filter1d(c_i_v, size=1), k=3)

        new_dif_mat_2[inum, :] = Spline(V1)

        Spline = UnivariateSpline(
            V0, uniform_filter1d(c_i_v, size=1), k=3)
        new_dif_mat_smooth_2[inum, :] = Spline(V1)

    spline_dif_mat_2 = new_dif_mat_2
    smooth_dif_mat = new_dif_mat_smooth_2

    'Plot the Results'
    

    if fold == 4:
    # if 1:
        t_axis = np.linspace(0, t_eval[-1], el_time, endpoint=True)  # cm -> um
        z_axis = np.linspace(0, Lp, el_num, endpoint=True)  # cm -> um

        Labels_X1 = 'time [sec]'  # mu: \u03BC
        Labels_X2 = 'z-direction [\u03BCm]'  # mu: \u03BC
        Labels_Y1 = 'Estimated Concentration [mM]'
        Labels_Y2 = 'Estimated DNA duplexes/UCN'
        Title = file_name

        fig = plt.figure()
        a0 = fig.add_axes([0, 0, 1, 1])
        a1 = a0.twinx()
        lines = []
        annotations = []
        print('spline_dif_mat=', np.shape(spline_dif_mat))
        print('smooth_dif_mat=',np.shape(smooth_dif_mat))

        for zi in z_position2:
            lines.append(a0.plot(t_axis, spline_dif_mat[:, zi]*10**6, '-',
                                 label="{0:.1f} μm".format(z_axis[zi]*10**4), linewidth=3, markersize=12)[0])
            annotations.append("{0:.1f} μm".format(z_axis[zi]*10**4))

            a1.plot(t_axis, spline_dif_mat[:, zi]*conc_conversion, '-',
                    linewidth=0.01, markersize=12)
            a1.plot(t_axis, smooth_dif_mat[:, zi]*conc_conversion, 'k--',
                    linewidth=3, markersize=12)
        idx = 12
        for line, label in zip(lines, annotations):
            x = t_axis
            # Index of closest x-value to x_pos
            idx -= 1
            # Get y-value of the line at the closest x-value
            y_pos = line.get_ydata()[idx]
            a0.annotate(label, xy=(xlim2[1]-(xlim2[1]-x[idx])/2, y_pos), xytext=(xlim2[1]-(xlim2[1]-x[idx])/2, y_pos+0.05e-8), color=line.get_color(),
                        fontsize=25, ha='center')

        line = Line2D([0, 1], [0, 1], linestyle='--', color='k', linewidth=3)
        a0.set_xlabel(Labels_X1)
        a0.set_ylabel(Labels_Y1, color='k')
        a0.set_ylim(ylim2[0], ylim2[1]*10**6)
        # a0.set_xlim([xlim2[0], xlim2[1]])
        a0.tick_params(colors='k', direction="in", length=4, width=2)
        leg2 = plt.legend([line], ['Spline Fit'],
                          loc='upper right',  shadow=True,
                          mode=None, fancybox=True)
        a1.set_ylim(ylim2[0], ylim2[1]*conc_conversion)
        a1.set_ylabel(Labels_Y2, color='k')
        a1.tick_params('y', colors='k', direction="in", length=4, width=2)
        plt.show()
        fig.savefig(save_file_path+str(xtal)+'_002.svg', format='svg', dpi=400, bbox_inches='tight')

        Title = file_name
        fig = plt.figure()
        a0 = fig.add_axes([0, 0, 1, 1])
        a1 = a0.twinx()
        lines = []
        annotations = []
        for ti in [0, 3, 7, 11, 15]:
            lines.append(a0.plot(z_axis*10**4, spline_dif_mat[ti, :]*10**6, '-',
                                 label=f"t={t_eval[ti]} sec", linewidth=3, markersize=12)[0])
            annotations.append(f"{t_eval[ti]} sec")

            a1.plot(z_axis*10**4, spline_dif_mat[ti, :]*conc_conversion, '-',
                    linewidth=0.01, markersize=12)
            a1.plot(z_axis*10**4, smooth_dif_mat[ti, :]*conc_conversion, 'k--',
                    linewidth=3, markersize=12)

        x_pos = xlim2[1]
        for line, label in zip(lines, annotations):
            x = z_axis*10**4
            # Index of closest x-value to x_pos
            idx = np.abs(x - x_pos).argmin()
            # Get y-value of the line at the closest x-value
            y_pos = line.get_ydata()[idx]
            a0.annotate(label, xy=(xlim2[1]-(xlim2[1]-x[idx])/2, y_pos), xytext=(xlim2[1]-(xlim2[1]-x[idx])/2, y_pos+0.05e-8), color=line.get_color(),
                        fontsize=25, ha='center')

        line = Line2D([0, 1], [0, 1], linestyle='--', color='k', linewidth=3)
        a0.set_xlabel(Labels_X2)
        a0.set_ylabel(Labels_Y1, color='k')
        a0.set_ylim(ylim2[0], ylim2[1]*10**6)
        a0.set_xlim([xlim2[0], xlim2[1]])
        a0.tick_params(colors='k', direction="in", length=4, width=2)
        leg2 = plt.legend([line], ['Spline Fit'],
                          loc='upper right',  shadow=True,
                          mode=None, fancybox=True)
        a1.set_ylabel(Labels_Y2, color='k')
        a1.set_ylim(ylim2[0], ylim2[1]*conc_conversion)
        a1.tick_params('y', colors='k', direction="in", length=4, width=2)
        plt.show()
        fig.savefig(save_file_path+str(xtal)+'_004.svg', format='svg', dpi=400, bbox_inches='tight')

    'Save data as .csv'
    z_axis = np.linspace(0, Lp, el_num, endpoint=True)  # cm -> um
    'Images data'
    csvfilename_1 = save_path_1 + '_fold' + \
        str(fold) + '_images_smooth_data.csv'

    csvfilename_2 = save_path_2 + '_fold' + \
        str(fold) + '_images_data.csv'
        
    with open(csvfilename_1, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(z_axis.reshape([1, el_num]))
        csvwriter.writerows(smooth_dif_mat)

    with open(csvfilename_2, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(z_axis.reshape([1, el_num]))
        csvwriter.writerows(spline_dif_mat)

    'Images BC data'
    csvfilename_1_r = save_path_1_r + '_fold_right' + \
        str(fold) + '_images_smooth_BC_data.csv'
    csvfilename_1_l = save_path_1_l + '_fold_left' + \
        str(fold) + '_images_smooth_BC_data.csv'

    csvfilename_2_r = save_path_2_r + '_fold_right' + \
        str(fold) + '_images_BC_data.csv'
    csvfilename_2_l = save_path_2_l + '_fold_left' + \
        str(fold) + '_images_BC_data.csv'

    with open(csvfilename_1_r, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(
            z_axis[:int(el_num/2)].reshape([1, int(el_num/2)]))
        csvwriter.writerows(smooth_dif_mat[:, int(el_num/2):][::-1])

    with open(csvfilename_1_l, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(
            z_axis[int(el_num/2):].reshape([1, int(el_num/2)]))
        csvwriter.writerows(smooth_dif_mat[:, :int(el_num/2)])

    with open(csvfilename_2_r, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(
            z_axis[:int(el_num/2)].reshape([1, int(el_num/2)]))
        csvwriter.writerows(spline_dif_mat[:, int(el_num/2):][::-1])

    with open(csvfilename_2_l, 'w') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(
            z_axis[:int(el_num/2)].reshape([1, int(el_num/2)]))
        csvwriter.writerows(spline_dif_mat[:, :int(el_num/2)])

print('Image detection completed.')
# %%
