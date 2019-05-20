# learning-ds-tutorial
This package includes demo scripts and a GUI simulation for learning stable non-linear Dynamical Systems (DS) from demonstrations using SEDS [1] and LPV-DS [2,3] approaches developed in LASA-EPFL.

### Running the demo scripts
There are three important demo scripts:
 - ```demo_SEDS.m```: Includes commented code blocks describing all the steps needed to learn a SEDS model on self-drawn data with a GUI or data loaded from the LASA handwriting dataset [1].
 - ```demo_LPVDS.m```: Includes commented code blocks describing all the steps needed to learn an LPV-DS model on self-drawn data with a GUI or data loaded from the LASA handwriting dataset [1] with different GMM fitting options and DS parameters optimization variants. 
- ```gui_learningDS.m```: Brings up a GUI including all options in the above script + a robot simulation (as shown below) . 
*A guided video explaining how to use the GUI, can be found in this [link](https://www.youtube.com/watch?v=5fQO9Oluih0)*

[![](https://github.com/nbfigueroa/learning-ds-tutorial/blob/master/img/GUI_2.png)](https://www.youtube.com/watch?v=5fQO9Oluih0)


**References**     
[1] Khansari Zadeh, S. M. and Billard, A. (2011) Learning Stable Non-Linear Dynamical Systems with Gaussian Mixture Models. IEEE Transaction on Robotics, vol. 27, num 5, p. 943-957.    
[2] Mirrazavi Salehian, S. S. (2018) Compliant control of Uni/ Multi- robotic arms with dynamical systems. PhD Thesis.  
[3] Figueroa, N. and Billard, A. (2018) A Physically-Consistent Bayesian Non-Parametric Mixture Model for Dynamical System Learning. In Proceedings of the 2nd Conference on Robot Learning (CoRL).

**Contact**: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387) (nadia.figueroafernandez AT epfl dot ch)
