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

' ======================================================================= '
' --- User Defined Settings --------------------------------------------- '
' ==========================/============================================= '

' - - STRING VARIABLES  - - - - - - - - - - - - - - - - - - - - - - - - - '
file_path = 'Crystals_images/'
file_name = 'Crystal1'    
xtal = 'C1'                         
# file_name = 'Crystal2'
# xtal = 'C2'
# file_name = 'Crystal3'
# xtal = 'C3'

fnameTiff = file_path + file_name + '.tiff'   # Multidimensional Tiff data file
fnameRGB = file_path + file_name + '.tif'     # RBG Tiff file for display
' - - NUMERIC VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - '
el_num = 100                # Constrain fit to this number of elements
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
    xlim1 = [-1, 67]     # x-axis range for fig0 (pixel)
    xlim2 = [-1, 47]     # x-axis range for fig1 (pixel)
    image0 = file_path+'Crystal1_0000.jpg'  # image path for plotting
    image1 = file_path+'Crystal1_0008.jpg'  # image path for plotting
    image2 = file_path+'Crystal1_0015.jpg'  # image path for plotting
    t_num = [0, 8, 15]   # Time index for image0, image1, and image 2
    trim_right = 152     # Pixel of right side of image to trim out for plotting
    trim_left = 152      # Pixel of left side of image to trim out for plotting
    trim_top = 100       # Pixel of top side of image to trim out for plotting
    trim_bottom = 120    # Pixel of bottom side of image to trim out for plotting

elif xtal == 'C2':
    trimmed_pixel = 40  # Pixel of image to trim out boundary optical effect
    height = 100        # Height of data slice (pixel)
    t_eval = [7, 127, 247, 367, 487, 607, 727, 847, 967,
              1087, 1207, 1327, 1447, 1567, 1687, 1807]  # Evaluation time points
    left = 140          # Start point of left side of the crystal (pixel)
    right = 250         # Start point of right side of the crystal (pixel)
    center = 210        # Center of the crystal (pixel)
    xlim1 = [-2, 53]    # x-axis range for fig0 (pixel)
    xlim2 = [-1, 36]    # x-axis range for fig1 (pixel)
    image0 = file_path+'Crystal2_0000.jpg'  # image path for plotting
    image1 = file_path+'Crystal2_0008.jpg'  # image path for plotting
    image2 = file_path+'Crystal2_0015.jpg'  # image path for plotting
    t_num = [0, 8, 15]  # Time index for image0, image1, and image 2
    trim_right = 123    # Pixel of right side of image to trim out for plotting
    trim_left = 123     # Pixel of left side of image to trim out for plotting
    trim_top = 90       # Pixel of top side of image to trim out for plotting
    trim_bottom = 70    # Pixel of bottom side of image to trim out for plotting

elif xtal == 'C3':
    trimmed_pixel = 15  # Pixel of image to trim out boundary optical effect
    height = 45         # Height of data slice             
    t_eval = [17, 137, 257, 377, 497, 617, 737, 857, 977,
              1097, 1217, 1337, 1457, 1577, 1697, 1817]    # Evaluation time points
    left = 83           # Start point of left side of the crystal (pixel)
    right = 204         # Start point of right side of the crystal (pixel)
    center = 190        # Center of the crystal (pixel)
    xlim1 = [-3, 155]   # x-axis range for fig0 (pixel)
    xlim2 = [-1, 110]   # x-axis range for fig1 (pixel)
    image0 = file_path+'Crystal3_0000.jpg' # image path for plotting
    image1 = file_path+'Crystal3_0008.jpg' # image path for plotting
    image2 = file_path+'Crystal3_0015.jpg' # image path for plotting
    t_num = [0, 8, 15]  # Time index for image0, image1, and image 2
    trim_right = 70     # Pixel of right side of image to trim out for plotting
    trim_left = 70      # Pixel of left side of image to trim out for plotting
    trim_top = 60       # Pixel of top side of image to trim out for plotting
    trim_bottom = 50    # Pixel of bottom side of image to trim out for plotting

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

    if fold == 4:
        fig0 = plt.figure()
        a0 = fig0.add_axes([0, 0, 1, 1])
        JDS_pt = [0, 3, 7, 11, 15]
        lines = []
        annotations = []
        for i in JDS_pt:
            lines.append(a0.plot(xtal_I_vec[i, :], '-',
                                 label=f"t={t_eval[i]} sec", linewidth=3, markersize=12)[0])
            annotations.append(f"{t_eval[i]} sec")

        x_pos = xlim1[1]
        for line, label in zip(lines, annotations):
            x = np.arange(0, np.shape(xtal_I_vec[i, :])[0], 1)
            # Index of closest x-value to x_pos
            idx = np.abs(x - x_pos).argmin()
            # Get y-value of the line at the closest x-value
            y_pos = line.get_ydata()[idx]
            a0.annotate(label, xy=(xlim1[1]-(xlim1[1]-x[idx])/2, y_pos), xytext=(xlim1[1]-(xlim1[1]-x[idx])/2, y_pos), color=line.get_color(),
                        fontsize=25, ha='center')

        plt.xlabel('pixel #')
        plt.ylabel('Intensity')
        plt.ylim([7000, 55000])
        plt.xlim([xlim1[0], xlim1[1]])
        plt.tick_params(direction="in", length=4, width=2)
        plt.show()


    'Convert from intensity values to concentrations based on standard'
    xtal_I = xtal_I_vec-I0
    xtal_I_vec_conc = 3.63042E-17*xtal_I**2 + 1.24742E-13*xtal_I + 3.60000E-10

    ' ======================================================================= '
    ' --- Fit Image Data to Diffusion Coef. --------------------------------- '
    ' ======================================================================= '
    Lp = round((xtal_r-xtal_l))*1e-4*pixel_length

    dif_mat = xtal_I_vec_conc
    Final_Dist = dif_mat[-1, :]
    x_dist = np.shape(dif_mat)[1]

    ' Matrix of concentrations in time vs position. '
    new_dif_mat = np.zeros([num_images, el_num])
    new_dif_mat_smooth = np.zeros([num_images, el_num])
    V0 = np.arange(1, x_dist+1) 
    V1 = np.linspace(1, x_dist, el_num, endpoint=True)

    for inum in range(0, num_images):
        c_i_v = dif_mat[inum, :]
        ' Interpolate to find values of c_i_v at coordinates V1.'
        # B-spline degree. Default is cubic, k=3.
        Spline = make_interp_spline(
            V0, uniform_filter1d(c_i_v, size=1), k=3)

        new_dif_mat[inum, :] = Spline(V1)

        Spline = UnivariateSpline(
            V0, uniform_filter1d(c_i_v, size=1), k=3)
        new_dif_mat_smooth[inum, :] = Spline(V1)

    spline_dif_mat = new_dif_mat
    smooth_dif_mat = new_dif_mat_smooth

    'Plot the Results'
    x_axis = np.linspace(0, Lp, el_num, endpoint=True)  # cm -> um

    if fold == 4:
        Labels_X = 'z-direction [\u03BCm]'  # mu: \u03BC
        Labels_Y1 = 'Concentration [mole/mL]'
        Labels_Y2 = 'DNA duplexes/UCN'
        Title = file_name

        fig0 = plt.figure()
        a0 = fig0.add_axes([0, 0, 1, 1])
        a1 = a0.twinx()
        lines = []
        annotations = []
        for i in JDS_pt:
            lines.append(a0.plot(x_axis*10**4, spline_dif_mat[i, :], '-',
                                 label=f"t={t_eval[i]} sec", linewidth=3, markersize=12)[0])
            annotations.append(f"{t_eval[i]} sec")

            a1.plot(x_axis*10**4, spline_dif_mat[i, :]*conc_conversion, '-',
                    linewidth=0.01, markersize=12)
            a1.plot(x_axis*10**4, smooth_dif_mat[i, :]*conc_conversion, 'k--',
                    linewidth=3, markersize=12)

        x_pos = xlim2[1]
        for line, label in zip(lines, annotations):
            x = x_axis*10**4
            # Index of closest x-value to x_pos
            idx = np.abs(x - x_pos).argmin()
            # Get y-value of the line at the closest x-value
            y_pos = line.get_ydata()[idx]
            a0.annotate(label, xy=(xlim2[1]-(xlim2[1]-x[idx])/2, y_pos), xytext=(xlim2[1]-(xlim2[1]-x[idx])/2, y_pos+0.05e-8), color=line.get_color(),
                        fontsize=25, ha='center')

        line = Line2D([0, 1], [0, 1], linestyle='--', color='k', linewidth=3)
        a0.set_xlabel(Labels_X)
        a0.set_ylabel(Labels_Y1, color='k')
        a0.set_ylim(0, 1e-7)
        a0.set_xlim([xlim2[0], xlim2[1]])
        a0.tick_params(colors='k', direction="in", length=4, width=2)
        leg2 = plt.legend([line], ['Spline Fit'],
                          loc='upper right',  shadow=True,
                          mode=None, fancybox=True)
        a1.set_ylabel(Labels_Y2, color='k')
        a1.set_ylim(0, 1e-7*conc_conversion)
        a1.tick_params('y', colors='k', direction="in", length=4, width=2)
        plt.show()

print('Image detection completed.')
# %%
