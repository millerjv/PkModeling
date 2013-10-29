# PkModeling
PkModeling is a [3D Slicer Version 4](http://www.slicer.org) Extension that provides pharmacokinetic modeling for dynamic contrast enhanced MRI (DCE MRI)[1][2].

PkModeling accepts volumetric timecourse data of signal intensities and computes parametric maps using either a two or three parameter Tofts model. The estimated parameters include

* Ktrans - volume transfer contrast constant between plasma and the extracellular-extravascular space at each voxel
* Ve - fractional volume for extracellular space at each voxel
* fpv - fractional plasma volume at each voxel
* MaxSlope - maximum slope of the of the time sequence at each voxel
* AUC - area under the curve of each voxel, measured from the bolus arrival time to the end of the time interval, normalized by the AUC of the arterial input function (AIF)
* R^2 - goodness of fit value. Since the parametric model is non-linear, R^2 is not strictly bounded by [-1,1]. But larger values still correspond to better fits.

PkModeling can also output a concentration curve view of the original volumetric timecourse as well as the "fitted" concentration curves resulting from the parametric model.

Estimation of the parametric model is controlled through a series of inputs including

* T1 Blood Value
* T1 Tissue Value - defaults to published value for prostate in de Bazelaire et al.[3]
* Relaxivity Value - default 0.0039 corresponds to the Gd-DPTA (Magnevist) at 3T, see Pintaske et al.[4]. This value needs to be adjusted for magnet strength and contrast agent.
* Hematocrit Value - volume percentage of red blood cells in blood
* AUC Time Interval Value - time interval for AUC calculation

Furthermore, an arterial input function (AIF) must be specified either by designating a mask corresponding to the voxels on which to base a patient specific estimate of the AIF, or by specifying a population derived AIF curve directly.

Finally, the estimation of the parametric maps can be restricted to a specified mask defining a region of interest.

Acquisition parameters relevent to the parametric model fitting are embedded in the input volumetric timecourse data, either as attributes on a NRRD file or extracted directly from the underlying DICOM structures

* TR Value - repetition time (in milliseconds)
* TE Value - echo time (in milliseconds)
* FA Value - flip angle (in degrees)
* Timestamps for the timecourses (in milliseconds)

# Visualization
See the [MultiVolumeExplorer](ttps://github.com/fedorov/MultiVolumeExplorer) module in the 3D Slicer.

# Authors
@millerjv, @fedorov, @zhuy

# References
[1]: Knopp MV, Giesel FL, Marcos H et al. "Dynamic contrast-enhanced magnetic resonance imaging in oncology." Top Magn Reson Imaging, 2001; 12:301-308.

[2]: Rijpkema M, Kaanders JHAM, Joosten FBM et al. "Method for quantitative mapping of dynamic MRI contrast agent uptake in human tumors." J Magn Reson Imaging 2001; 14:457-463.

[3]: de Bazelaire, C.M., et al. "MR imaging relaxation times of abdominal and pelvic tissues measured in vivo at 3.0 T: preliminary results." Radiology, 2004. 230(3): p. 652-9.

[4]: Pintaske J, Martirosian P, Graf H, Erb G, Lodemann K-P, Claussen CD, Schick F. "Relaxivity of Gadopentetate Dimeglumine (Magnevist), Gadobutrol (Gadovist), and Gadobenate Dimeglumine (MultiHance) in human blood plasma at 0.2, 1.5, and 3 Tesla." Investigative radiology. 2006 March;41(3):213–21.

