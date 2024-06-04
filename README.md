Biophysics
==========
This repository incorporates several tools that enable the analysis of transportation of 8mer dsDNA inside CJ protein porous crystal from a time-lapse confocal images.

## Scripts Description

### 1. Image_2D_smoothing_20_segments.py

Process time-lapse confocal images.

Starts from defining the crystal edges in the image that is stored in Crystals_images folder, followed by segemented data into 20 sections across z- and y-direction of tranportation. Finally, converted light intensity to estimated dsDNA concentration from 3-parameters standard curve equation and save data into .csv file and is stored in Crystal_csv and Crystals_smooth_csv folder.

### 2. Fick_Law_Implementation.py

Calculating a range of diffusion rate constant through Fick's law.

Starts from calculating accumulation rates (RA) and driving force (DF) for all position and all time points for each segmented data from folder Crystals_smooth_csv, and stored into folder Transport_tables. Followed by plotting RA vs DF and calculate the slope as diffusion rate constant. Finally, store all diffusion rate constant along with its position and time points to .csv file in folder Whisker_plot_csv.

### 3. Whisker_plot.py

Alanlysis of diffusion rate constants from Fick's law implementation.

Plot data of diffusion rate constant from folder Whisker_plot_csv as whisker plot and calculate the medium value of diffusion rate constant for each segmented data.

### 4. FVM.py

Validate the medium diffusion rate constant by fitting data with finite volume model (FVM). Included two files, FVM_1D_Dl_Dr.py and FVM_Function.py.

Assigned medium diffusion rate values to parameteres of Dl (Diffusion rate of left edge) and Dr (Diffusion rate of left edge). Run the script FVM_1D_Dl_Dr.py by using finite volume model (FVM) and boundary concentration (folder Crystal_csv) as a function of time. Finally, calculate sum of squared deviation (SSD) of each segmented data fitting to FVM.
