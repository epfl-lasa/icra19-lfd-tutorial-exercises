To start the demo:

1- Execute the script setup_demo.m

2- Run modulated_ds_interface.m

   Usage: 1: modulated_ds_interface
          2: modulated_ds_interface(generatePerturbation, makeAMovie, loadSavedSurfaceModel)
          generatePerturbation: (Boolean) Generate a perturbation on the robot during the
                                          simulation. Default = false.
          makeAMove: (Boolean) Generate a movie of the simulation. Default = false.
          loadSavedSurfaceModel: (Boolean) Load the surface model previously saved.
                                 Default = false.



Nadia Comments:
- Current code is not stand-alone. Depends on svm_regressor function from ML_toolbox [Can be added]
- Current code will only work for MATLAB version 2017a upwards, due to the vecnorm() function, I have a fix for this in my lags code
- gpml folder not necessary [Removed]
- Some of the functions for the GUI depend on Image Processing Toolbox, this might be a problem for people who don't have that toolbox installed
