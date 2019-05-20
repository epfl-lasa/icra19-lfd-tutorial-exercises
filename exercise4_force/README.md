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
- Current code is not stand-alone. Depends on svm_regressor function from ML_toolbox
- Current code will only work for MATLAB version 2017a upwards, due to the vecnorm() function
- gpml folder not necessary [Removed]