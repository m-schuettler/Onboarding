## for running 'decode_barcode()' function on one fov, to be run by 'run_run.sh'

import glob
import os
import functions as fct # import functions.py as module
import sys              # import functions needed to run in 'run_run.sh'

# define 'fov' list
fov = []
searchpattern = '../data/selected-tiles/selected-tiles/out_opt_flow_registered*c01_DAPI.tif'  # get all names of fovs
for file in glob.glob(searchpattern):
    fovname = os.path.basename(file).removeprefix('out_opt_flow_registered').removesuffix('c01_DAPI.tif')
    fov.append(fovname)
b = int(sys.argv[1])        # 'b' is defined when running 'bashdebug.sh'
fov = fov[b]                # specify which fov should be used right now
print('Current fov: ' + str(b+1), end='\r')

subfolder = sys.argv[2]     # 'subfolder' is defined when running 'run_run.sh'

# run 'decode_barcode()' function
fct.decode_barcode(fov, subfolder)