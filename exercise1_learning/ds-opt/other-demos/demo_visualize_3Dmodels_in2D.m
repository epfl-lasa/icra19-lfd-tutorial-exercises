%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Load one of the 3D Datasets  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc
%%%%%%%%%%%%%%%%%% Select an Axis-Aligned Dataset %%%%%%%%%%%%%%%%%%%%%%%%%%
% 6:  Via-point Dataset     (3D) * 9  trajectories recorded at 100Hz
% 7:  Sink Dataset          (3D) * 11 trajectories recorded at 100Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg_dir         = '/home/nbfigueroa/Dropbox/PhD_papers/CoRL-2018/code/ds-opt/';
chosen_dataset  = 7; 
sub_sample      = 2; % '>2' for real 3D Datasets, '1' for 2D toy datasets
nb_trajectories = 9; % For real 3D data only
[Data, Data_sh, att, x0_all, data, dt] = load_dataset_DS(pkg_dir, chosen_dataset, sub_sample, nb_trajectories);

% Position/Velocity Trajectories
vel_samples = 10; vel_size = 0.5; 
[h_data, h_att, h_vel] = plot_reference_trajectories_DS(Data, att, vel_samples, vel_size);

% Extract Position and Velocities
M          = size(Data,1)/2;    
Xi_ref     = Data(1:M,:);
Xi_dot_ref = Data(M+1:end,:);   

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Load pre-learned lpv-DS model from Mat file  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DS_name = '/3D-Sink/3D-Sink_pqlf_2';
% DS_name = '3D-Via-point_pqlf_2';
matfile = strcat(pkg_dir,'/models/', DS_name,'.mat');
load(matfile)
if constr_type == 1
    ds_lpv = @(x) lpv_ds(x-repmat(att,[1 size(x,2)]), ds_gmm, A_g, b_g);
else
    ds_lpv = @(x) lpv_ds(x, ds_gmm, A_k, b_k);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Visualize DS Streamline in 2D axis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot YZ coordinates of DS
yz = [2 3];
ds_plot_options = [];
ds_plot_options.sim_traj   = 1;      % To simulate trajectories from x0_all
ds_plot_options.x0_all     = x0_all; % Initial Points for reproductions
ds_plot_options.dimensions = yz;
ds_plot_options.attractor  = att;
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_lpv, ds_plot_options);
title('YZ-slice of learned DS', 'Interpreter','LaTex','FontSize',20)

%% Plot YX coordinates of DS
yx = [2 1];
ds_plot_options = [];
ds_plot_options.sim_traj   = 1;      % To simulate trajectories from x0_all
ds_plot_options.x0_all     = x0_all; % Initial Points for reproductions
ds_plot_options.dimensions = yx;
ds_plot_options.attractor  = att + [0 0 0.005]';
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_lpv, ds_plot_options);
title('YX-slice of learned DS', 'Interpreter','LaTex','FontSize',20)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Load pre-learned lpv-DS model from Mat file  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DS_name = '/3D-Sink/3D-Sink_seds';
matfile = strcat(pkg_dir,'/models/', DS_name,'.mat');
load(matfile)
ds_seds = @(x) GMR_SEDS(Priors,Mu,Sigma,x-repmat(att,[1 size(x,2)]),1:M,M+1:2*M);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Visualize DS Streamline in 2D axis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot YZ coordinates of DS
yz = [2 3];
ds_plot_options = [];
ds_plot_options.sim_traj   = 1;      % To simulate trajectories from x0_all
ds_plot_options.x0_all     = x0_all; % Initial Points for reproductions
ds_plot_options.dimensions = yz;
ds_plot_options.attractor  = att;
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_seds, ds_plot_options);
title('YZ-slice of learned DS', 'Interpreter','LaTex','FontSize',20)

%% Plot YX coordinates of DS
yx = [1 2];
ds_plot_options = [];
ds_plot_options.sim_traj   = 1;      % To simulate trajectories from x0_all
ds_plot_options.x0_all     = x0_all; % Initial Points for reproductions
ds_plot_options.dimensions = yx;
ds_plot_options.attractor  = att + [0 0 0.05]';
[hd, hs, hr, x_sim] = visualizeEstimatedDS(Xi_ref, ds_seds, ds_plot_options);
title('XY-slice of learned DS', 'Interpreter','LaTex','FontSize',20)
