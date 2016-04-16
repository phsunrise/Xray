# Xray
## Data analysis procedure:
1. Create "parameters.pickle" in 'w' (not 'wb') mode, containing the data directory: {'data_dir': 'path/to/data/directory'}
2. Run "locate_pads.py" to locate the subpads. May need to look for a good run and change contrast factor A
3. Run "calibration.py" to do the calibration. May need to change the twotheta_deg array inside "calibraiton.py" to match the diffraction pattern. The array needs to be in increasing order of 2theta
4. Run "process.py" for data processing.

##Command line arguments:
* -p: pad number
* -r: run number
* -i: image number
* -b: background run number
* -d: debug mode (for "process.py", will show subpad offset plots and exit before fitting)
