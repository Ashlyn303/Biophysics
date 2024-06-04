# %%
from combine_svgs_v3 import combine_svgs
from glob import glob
import os

import numpy as np
output_path = os.getcwd()
print('output_path=', output_path)

dirs = glob('0*')
dirs.sort()
print('There are {} dir(s).'.format(len(dirs)))
dir_num = 0
for d in dirs:
    dir_num += 1
    print(str(dir_num)+') ' + d)

# Define matrix-like shape for svg group
shapes = [[6, 2], [9, 3], [6, 4], [1, 1], [6, 6]]
# Define figure order for each svg into svg group
# It allows to have blank space in between the svg group
# Also to add lables for each svg file
orders = [[1, 3, 2, 4, 5, 7, 6, 8, 9, 11, 10, 12],
          [1, 2, 4, 5, 6, 7, 8, 9, 
           10, 11, 13, 14, 15, 16, 17, 18,
           19, 20, 22, 23, 24, 25, 26, 27],
          range(24),
          [1],
          range(36)]

labels = [['B', 'C', 'D', 'E',
           'B', 'C', 'D', 'E',
           'B', 'C', 'D', 'E'],
          ['A', 'B', 'C Left Side', 'D Left Side', 'E Anomalous Map', 'F Right Side', 'G Right Side', 'H Anomalous Map',
           'A', 'B', 'C Left Side', 'D Left Side', 'E Anomalous Map', 'F Right Side', 'G Right Side', 'H Anomalous Map',
           'A', 'B', 'C Left Side', 'D Left Side', 'E Anomalous Map', 'F Right Side', 'G Right Side', 'H Anomalous Map'],
          ['A', '', '', '', 'B', '', '', '',
           'A', '', '', '', 'B', '', '', '',
           'A', '', '', '', 'B', '', '', ''],
          [''],
          ['A', '', '', '', '', '', 'B', '', '', '', '', '',
           'A', '', '', '', '', '', 'B', '', '', '', '', '',
           'A', '', '', '', '', '', 'B', '', '', '', '', '']]

# iterate through each subfolder and create svg group file
dir_num = 0
for d in dirs:
    dir_num += 1
    os.system('cp combine_svgs_v3.py %s/' % d)
    os.chdir(d)

    svg_paths = glob('*.svg')
    svg_paths.sort()
    print('There are {} svg(s) in {}.'.format(len(svg_paths), d))
    _ = 0
    for svg in svg_paths:
        _ += 1
        print(str(_)+') ' + svg)

    shape = shapes[dir_num-1]
    order = orders[dir_num-1]
    label = labels[dir_num-1]
    output_path = '/Users/szu-hsuan/Dropbox/python3/Biophysics_paper/savefig' + \
        '/combined_svg_00{}.svg'.format(dir_num)
    combine_svgs(svg_paths, shape, order, label, output_path)
    os.chdir('../')

# %%
