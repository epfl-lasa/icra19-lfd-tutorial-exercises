# exercise1_learning
This package includes GUI simulations and demo scripts for learning stable non-linear Dynamical Systems (DS) from demonstrations using LPV-DS [1] and SEDS [2] approaches developed in LASA-EPFL.

### Running the gui for exercise sessions
- ```gui_lpvDS.m```: Brings up a GUI including for estimation of Globally Asymptotically Stable Dynamical Systems focused on Linear Parameter Varying formulation with GMM-based mixing function and different Lyapunov candidate functions as proposed in [1], where a non-linear DS formulated as:
<p align="center">
<img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/img/f_x.gif"></>
                                                                           
is learned from demonstrations in a decoupled manner. Where the GMM parameters <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/img/theta_gamma.gif"> used to parametrize <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/img/gamma.gif"> are estimated via the physically-consistent GMM approach proposed in [1] and provided in [phys-gmm](https://github.com/nbfigueroa/phys-gmm) and the DS parameters <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/img/DS_params.gif"> are estimated via semi-definite optimization problem that ensures global asymptotic stability of the system via constraints derived from either a:
- QLF (Quadratic Lyapunov Function): <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/img/stab_qlf.gif">
- P-QLF(Parametrized QLF):  <img src="https://github.com/nbfigueroa/ds-opt/blob/master/figs/img/stab_pqlf.gif">
[![](https://github.com/nbfigueroa/learning-ds-tutorial/blob/master/img/GUI_2.png)](https://www.youtube.com/watch?v=5fQO9Oluih0)
*A guided video explaining how to use the GUI, can be found in this [link](https://www.youtube.com/watch?v=5fQO9Oluih0)*


- [Optional] ```gui_seDS.m```: Brings up a GUI including SEDS (constrained-GMR) learning approach with robot simulation. 
[![](https://github.com/nbfigueroa/learning-ds-tutorial/blob/master/img/GUI_2.png)](https://www.youtube.com/watch?v=5fQO9Oluih0)
*A guided video explaining how to use the GUI, can be found in this [link](https://www.youtube.com/watch?v=5fQO9Oluih0)*


### Using algorithms on own data
Within the [ds_opt](https://github.com/epfl-lasa/icra19-lfd-tutorial-exercises/tree/master/exercise1_learning/ds-opt) folder, we provide demo scripts for using any of the presented learning algorithms on test 2D and 3D data. These scripts are useful if you want to learn DS with your own datasets. 


**References**     
[1] Figueroa, N. and Billard, A. (2018) A Physically-Consistent Bayesian Non-Parametric Mixture Model for Dynamical System Learning. In Proceedings of the 2nd Conference on Robot Learning (CoRL).
[2] Khansari Zadeh, S. M. and Billard, A. (2011) Learning Stable Non-Linear Dynamical Systems with Gaussian Mixture Models. IEEE Transaction on Robotics, vol. 27, num 5, p. 943-957.    


**Contact**: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387) (nadia.figueroafernandez AT epfl dot ch)
