Biophysics
==========
This repository incorporates several tools that enable the analysis of transportation of 8mer dsDNA inside CJ protein porous crystal from a time-lapse confocal images.
# Scripts Description
**1. Image_detect_20_segments.py**

Process time-lapse confocal images.

Starts from defining the crystal edges in the image, followed by segemented data into 20 sections across z- and y-direction of tranportation. Finally, converted light intensity to estimated dsDNA concentration from 3-parameters standard curve equation and save data into .csv file.
**2. Fick_Law_Implementation.py**

Calculating a range of diffusion rate constant through Fick's law.

Starts from calculating accumulation rates (RA) and driving force (DF) for all position and all time points for each segmented data. Followed by plotting RA vs DF and calculate the slope as diffusion rate constant. Finally, store all diffusion rate constant along with its position and time points to .csv file.
**3. Whisker_plot.py**

Alanlysis of diffusion rate constants from Fick's law implementation.

Plot data os diffusion rate constant as whisker plot and calculate the medium value of diffusion rate constant for each segmented data.
**4. FVM.py**

Validate the medium diffusion rate constant by fitting data with finite volume model (FVM). Included two files, FVM_1D_Dl_Dr.py and FVM_Function.py.

