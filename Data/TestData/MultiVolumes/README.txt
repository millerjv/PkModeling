The nrrd files in this folder is generated to test PkModeling, they are all 4D images.
Size: 10*10*5*40

Parameters of synthetic test data for testing SignalIntensitiesToConcentrationValues
FA = 15 degree
TR = 3.984
T1 = 1600
s0 = 10 (same as in the input data)
relaxivity = 0.0049

SyntheticMultiVolumeInput.nrrd: signal intensity values, each voxel has the following values in time domain
10, 10, 10, 10, 10, 10, 10, 10, 50, 100, 50, 10, 10, 10, 10,10, 10, 10, 10, 10,10, 10, 10, 10, 10,10, 10, 10, 10, 10,10, 10, 10, 10, 10,10, 10, 10, 10, 10,10, 10, 10, 10, 10,10, 10, 10, 10, 10,
SyntheticMultiVolumeOutput.nrrd: concentration values calculated from the signal intensity values, each voxel has the following values in time domain
0, 0, 0, 0, 0, 0, 0, 0, 0.7673, 3.4819, 0.7673, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
SyntheticMultiVolumeAllZero.nrrd: all zero values
SyntheticMultiVolumeAllTen.nrrd: all "10" values
SyntheticMultiVolumeAllNagTen.nrrd: all "-10" values
