%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Code for incremental LPV-DS Learning Algorithm from paper:         %
%  'A Physically-Consistent Bayesian Non-Parametric Mixture Model for     %
%   Dynamical System Learning.'; N. Figueroa and A. Billard; CoRL 2018    %
% With this script you can load 2D toy trajectories or even real-world    %
% trajectories acquired via kinesthetic taching and test the different    %
% GMM fitting approaches.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 Learning Algorithms and Systems Laboratory,          %
% EPFL, Switzerland                                                       %
% Author:  Nadia Figueroa                                                 % 
% email:   nadia.figueroafernandez@epfl.ch                                %
% website: http://lasa.epfl.ch                                            %
%                                                                         %
% This work was supported by the EU project Cogimon H2020-ICT-23-2014.    %
%                                                                         %
% Permission is granted to copy, distribute, and/or modify this program   %
% under the terms of the GNU General Public License, version 2 or any     %
% later version published by the Free Software Foundation.                %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General%
% Public License for more details                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 1: Draw with GUI %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear all Data
close all; clear all; clc;
%% Draw batches of data
fig1 = figure('Color',[1 1 1]);
limits = [-6 6 -2 0.5];
axis(limits)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.25, 0.55, 0.2646 0.4358]);
grid on

% Global Attractor of DS
att_g = [0 0]';
radius_fun = @(x)(1 - my_exp_loc_act(5, att_g, x));
scatter(att_g(1),att_g(2),100,[0 0 0],'d'); hold on;

% Draw Reference Trajectories
data = draw_mouse_data_on_DS(fig1, limits);
Data = []; x0_all = [];
for l=1:length(data)    
    % Check where demos end and shift
    data_ = data{l};
    data_(1:2,:) = data_(1:2,:) - repmat(data_(1:2,end), [1 length(data_)]);
    data_(3:4,end) = zeros(2,1);
    Data = [Data data_];
    x0_all = [x0_all data_(1:2,1)];
end

% Position/Velocity Trajectories
Xi_ref     = Data(1:2,:);
Xi_dot_ref = Data(3:end,:);

%% Store as first batch
Data_1 = Data;
x0_all_1 = x0_all;
Xi_ref_1 = Xi_ref;
Xi_dot_ref_1 = Xi_dot_ref;

%% Store as second batch
Data_2 = Data;
x0_all_2 = x0_all;
Xi_ref_2 = Xi_ref;
Xi_dot_ref_2 = Xi_dot_ref;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 2: Choose from LASA DATASET %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
% clear all; close all; clc

% Global Attractor of DS
att_g = [0 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load Chosen Dataset %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-Models are Motions 27-28-29-30...
m = 28;
sample = 2;
show_trajectories = 1;
[demos, name] = load_LASA_dataset_shape(m);

% Global Attractor of DS
att_g = [0 0]';

% Extracting data for model 1
Data_1 = []; x0_all_1 = [];
for l=5:7
    % Check where demos end and shift
    data_ = [demos{l}.pos(:,1:sample:end); demos{l}.vel(:,1:sample:end);];
    Data_1 = [Data_1 data_];
    x0_all_1 = [x0_all_1 data_(1:2,20)];
    clear data_
end

% Position/Velocity Trajectories
Xi_ref_1     = Data_1(1:2,:);
Xi_dot_ref_1 = Data_1(3:end,:);


% Extracting data for model 2
Data_2 = []; x0_all_2 = [];
for l=1:4
    % Check where demos end and shift
    data_ = [demos{l}.pos(:,1:sample:end); demos{l}.vel(:,1:sample:end);];
    Data_2 = [Data_2 data_];
    x0_all_2 = [x0_all_2 data_(1:2,20)];
    clear data_
end

% Position/Velocity Trajectories
Xi_ref_2     = Data_2(1:2,:);
Xi_dot_ref_2 = Data_2(3:end,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA LOADING OPTION 3: Load from Pre-drawn Examples %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example Figure. 6
clear all; close all; clc
load('./datasets/CoRL-paper-Datasets/incremental_1.mat')

% Example Figure. 7
clear all; close all; clc
load('./datasets/CoRL-paper-Datasets/incremental_2.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Show Trajectories either from Option 1 or 2       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_trajectories = 1;
if show_trajectories
    figure('Color',[1 1 1]);
    scatter(Xi_ref_1(1,:), Xi_ref_1(2,:),10,[1 0 0],'filled'); hold on
    scatter(Xi_ref_2(1,:), Xi_ref_2(2,:),10,[1 0 1],'filled'); hold on
    grid on; box on;
    xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
    ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
    legend({'Batch $b$', 'Batch $b+1$'},'Interpreter','LaTex','FontSize',15)
    title('Reference Trajectories','Interpreter','LaTex','FontSize',20);
    xLimits = get(gca,'XLim'); 
    yLimits = get(gca,'YLim');
    limits = [xLimits yLimits];
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 2: Learn LPV-DS for batch 1 %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% GMM Estimation Algorithm %%%%%%%%%%%%%%%%%%%
% 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
% 1: GMM-EM Model Selection via BIC
% 2: GMM via Competitive-EM
% 3: CRP-GMM via Collapsed Gibbs Sampler
gmm_type = 0;
%%%%%%%%%%%%%%%%%%% DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%
% Type of constraints/optimization 
constr_type = 0;      % 0:'convex':     A' + A < 0
                      % 1:'non-convex': A'P + PA < 0
                      % 2:'non-convex': A'P + PA < -Q given P                                  
init_cvx    = 0;      % 0/1: initialize non-cvx problem with cvx                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Learn model for batch 1
[ds_gmm_1, A_g_1, b_g_1, P_est_1] = learn_LPVDS(gmm_type,constr_type, init_cvx, Data_1, att_g, limits);

% Add number of points used to estimate this GMM
ds_gmm_1.N = size(Data_1,2);

% Create DS function handle
ds_lpv_1 = @(x) lpv_ds(x, ds_gmm_1, A_g_1, b_g_1);

%% Plot DS batch 1
fig1 = figure('Color',[1 1 1]);
[hd] = scatter(Xi_ref_1(1,:),Xi_ref_1(2,:),20,[1 0 0],'filled'); hold on
[hs] = plot_ds_model(fig1, ds_lpv_1, att_g, limits,'medium'); hold on;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
% legend({'Batch $b$'},'Interpreter','LaTex','FontSize',15)

% Simulate trajectories and plot them on top
plot_repr = 0;
if plot_repr
    opt_sim = [];
    opt_sim.dt = 0.01;
    opt_sim.i_max = 3000;
    opt_sim.tol = 0.1;
    opt_sim.plot = 0;
    [x_lpv xd_lpv]=Simulation(x0_all_1 ,[],ds_lpv_1, opt_sim);
    [hr] = scatter(x_lpv(1,:),x_lpv(2,:),5,[0 0 0],'filled'); hold on    
end
[~, est_labels_1] =  my_gmm_cluster(Xi_ref_1, ds_gmm_1.Priors, ds_gmm_1.Mu, ds_gmm_1.Sigma, 'hard', []);
plotGMMParameters( Xi_ref_1, est_labels_1, ds_gmm_1.Mu, ds_gmm_1.Sigma,fig1);
title('Initial batch $\theta_{\gamma}^b$ and $\mathbf{f}(\xi)^b$', 'Interpreter','LaTex','FontSize',20)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 3: Learn CP-GMM for batch 2 %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GMM Estimation Algorithm %%%%
% 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
% 1: GMM-EM Model Selection via BIC
% 2: GMM via Competitive-EM - rarely works..
% 3: CRP-GMM via Collapsed Gibbs Sampler

est_options = [];
est_options.type        = gmm_type; % GMM Estimation Alorithm Type    
est_options.maxK        = 15;       % Maximum Gaussians for Type 1/2
est_options.do_plots    = 1;        % Plot Estimation Statistics
est_options.adjusts_C   = 1;        % Adjust Sigmas
est_options.fixed_K     = 5;        % Fix K and estimate with EM
est_options.exp_scaling = 1;        % Scaling for the similarity to improve locality

% Discover Local Models
sample = 2;
[Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref_2(:,1:sample:end), Xi_dot_ref_2(:,1:sample:end), est_options);

%% Extract Cluster Labels
est_K_2      = length(Priors0); 
Priors_2 = Priors0; Mu_2 = Mu0; Sigma_2 = Sigma0;
[~, est_labels_2] =  my_gmm_cluster(Xi_ref_2, Priors_2, Mu_2, Sigma_2, 'hard', []);
clear ds_gmm_2; 
ds_gmm_2.Mu = Mu_2; ds_gmm_2.Sigma = Sigma_2; 
ds_gmm_2.Priors = Priors_2; 
% Adjust Covariance Matrices
if est_options.adjusts_C  == 1
    tot_scale_fact = 1; rel_scale_fact = 0.15;
    Sigma = adjust_Covariances(Sigma0, tot_scale_fact, rel_scale_fact);
    ds_gmm_2.Sigma = Sigma;
end  
% Add number of points used to estimate this GMM
ds_gmm_2.N = size(Data_2,2);

% Visualize Cluster Parameters on Manifold Data
plotGMMParameters( Xi_ref_2, est_labels_2, Mu_2, Sigma_2);
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
title('Batch $b+1$ $\mathcal{PC}-$GMM $\theta_{\gamma}^{b+1}$','Interpreter','LaTex','FontSize',20)
ml_plot_gmm_pdf(Xi_ref, Priors_2, Mu_2, Sigma_2, limits)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 4: Check for GMM Similarities of batch 2 wrt. batch 1 %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_K_1 = length(ds_gmm_1.Priors);
clc;
prune_components = [];
for k=1:est_K_2
   fprintf('********** Evaluating component k=%d from GMM bathc b+1 **********\n',k) 
   D_k = Xi_ref_2(:,est_labels_2==k);
   d = size(D_k,1);  n = size(D_k,2);  
   nu = d*(d+1)/2; 
   
   % Variables of j-th component
   Mu_k    = ds_gmm_2.Mu(:,k);
   Sigma_k = ds_gmm_2.Sigma(:,:,k);

   for j=1:est_K_1       
       % Variables of j-th component
       Mu_j    = ds_gmm_1.Mu(:,j); 
       Sigma_j = ds_gmm_1.Sigma(:,:,j);
       
       % Test KL-Divergence
       KL_div_1 = 0.5 * (log(det(Sigma_k)/det(Sigma_j)) - d + trace(inv(Sigma_k)*Sigma_j) + ...
                      (Mu_k-Mu_j)'*inv(Sigma_k)*(Mu_k-Mu_j));
       KL_div_2 = 0.5 * (log(det(Sigma_j)/det(Sigma_k)) - d + trace(inv(Sigma_j)*Sigma_k) + ...
                      (Mu_j-Mu_k)'*inv(Sigma_j)*(Mu_j-Mu_k));
       
       if KL_div_1 < 1 || KL_div_2 < 1
            fprintf(2,'KL-divergence(j,k)=%2.2f and KL-divergence(k,j)=%2.2f\n',KL_div_1, KL_div_2);
            fprintf(2,'Kb(%d) and Kb+1(%d) are practically the same! \n',j, k);
            prune_ = [k j];
            prune_components = [prune_components; prune_];
       else
            fprintf('KL-divergence(j,k)=%2.2f and KL-divergence(k,j)=%2.2f\n',KL_div_1, KL_div_2);
       end
   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 4: Prune-Merge/Concatenate GMMs and Estimate DS parameters %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(prune_components)
    fprintf(2,'*** No Pruning Necessary! Merging GMM components of both batches! ***\n');
    merged_ds_gmm = mergeGMMs(ds_gmm_1, ds_gmm_2);
    %%%%%%%%  LPV system sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%            
    [A_g_2, b_g_2, ~] = optimize_lpv_ds_from_data(Data_2, att_g, 0, ds_gmm_2, [], 0);          
    
else
    fprintf(2,'*** Need to prune components from batch b+1! ***\n');    
    % Before pruning.. estimate DS parameters
    %%%%%%%%  LPV system sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%            
    [A_g_2_, b_g_2_, ~] = optimize_lpv_ds_from_data(Data_2, att_g, 0, ds_gmm_2, [], 0);           
       
    remove_batch_2 = prune_components(:,1);
    update_batch_1 = prune_components(:,2);    
    K_prune = size(prune_components,1);

    % Auxiliary variable needed for computation
    M_k = ds_gmm_2.Priors*ds_gmm_2.N;
    M = ds_gmm_2.N;
    N = ds_gmm_1.N;
    
    K_batch1     = length(ds_gmm_1.Priors);
    K_batch2     = length(ds_gmm_2.Priors);    
    K_merge      = K_batch1  + K_batch2  - K_prune;    
    Mu_merge     = zeros(2,K_merge);
    Sigma_merge  = zeros(2,2,K_merge);
    Priors_merge = zeros(1,K_merge);
    
    % Fill in parameters for K_batch1
    for k=1:K_batch1
        if sum(k == update_batch_1)
            id_other_comp = remove_batch_2(update_batch_1 == k);
            fprintf('updating %d-th component of b with %d-th component b+1 \n',k, id_other_comp);
            Mu_j    = ds_gmm_1.Mu(:,k);
            Sigma_j = ds_gmm_1.Sigma(:,:,k);
            Priors_j = ds_gmm_1.Priors(k);
            
            Mu_k    = ds_gmm_2.Mu(:,id_other_comp);
            Sigma_k = ds_gmm_2.Sigma(:,:,id_other_comp);
            
            Mu_merge(:,k) = (N*Priors_j*Mu_j + M_k(id_other_comp)*Mu_k)/(N*Priors_j + M_k(id_other_comp));
            Sigma_merge(:,:,k) = (N*Priors_j*Sigma_j + M_k(id_other_comp)*Sigma_k)/(N*Priors_j + M_k(id_other_comp)) + ...
                                  (N*Priors_j*(Mu_j*Mu_j') + M_k(id_other_comp)*(Mu_k*Mu_k'))/(N*Priors_j + M_k(id_other_comp)) - (Mu_merge(:,k)*Mu_merge(:,k)');
            Priors_merge(k)  = (N*Priors_j + M_k(id_other_comp)) / (N + M);
            
        else
            fprintf('copying %d-th component\n',k)
            Mu_merge(:,k) = ds_gmm_1.Mu(:,k);
            Sigma_merge(:,:,k) = ds_gmm_1.Sigma(:,:,k);
            Priors_merge(k)  = (N*ds_gmm_1.Priors(k)) / (N + M);
        end
        
    end
        
    pruned_gmm = ds_gmm_2;
    idx_= 1:K_batch2; 
    pruned_idx    = setdiff(idx_,remove_batch_2');
    Mu_pruned     = ds_gmm_2.Mu(:,pruned_idx);
    Sigma_pruned  = ds_gmm_2.Sigma(:,:,pruned_idx);
    Priors_pruned = M_k(pruned_idx)/(N+M);

    % Add new components
    Mu_merge(:,K_batch1+1:end) = Mu_pruned;
    Sigma_merge(:,:,K_batch1+1:end) = Sigma_pruned;
    Priors_merge(K_batch1+1:end) = Priors_pruned;
    
    % Create new GMM structure
    merged_ds_gmm = [];
    merged_ds_gmm.Mu = Mu_merge; merged_ds_gmm.Sigma = Sigma_merge; merged_ds_gmm.Priors = Priors_merge;
    merged_ds_gmmN = N + M;
   
    % Extract desired A,b's
    A_g_2 = A_g_2_(:,:,pruned_idx);
    b_g_2 = b_g_2_(:,pruned_idx);
end

% Visualize Cluster Parameters on Manifold Data
[~, est_labels_merged] =  my_gmm_cluster([Xi_ref_1 Xi_ref_2], merged_ds_gmm.Priors, merged_ds_gmm.Mu, merged_ds_gmm.Sigma, 'hard', []);
plotGMMParameters( [Xi_ref_1 Xi_ref_2], est_labels_merged, merged_ds_gmm.Mu, merged_ds_gmm.Sigma);
title('Updated batch $b+1$ $\tilde{\theta}_{\gamma}^{b+1}$', 'Interpreter','LaTex','FontSize',20)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%  Step 5: Concatenate DS and Plot resulting DS %%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Concatenate f^b(\xi) and f^b+1(\xi) %%%%%%%%
K_merged = length(merged_ds_gmm.Priors);
A_merged = zeros(2,2,K_merged);
b_merged = zeros(2,K_merged);

% Initial Set of DS
A_merged(:,:,1:est_K_1)     = A_g_1;
A_merged(:,:,est_K_1+1:end) = A_g_2;
b_merged(:,1:est_K_1)     = b_g_1;
b_merged(:,est_K_1+1:end) = b_g_2;

% Create DS function handle
ds_lpv_merged = @(x) lpv_ds(x, new_ds_gmm, A_new, b_new);

%% %%%%%%%%%%%%    Plot Resulting Merged DS  %%%%%%%%%%%%%%%%%%
fig1 = figure('Color',[1 1 1]);
[hd] = scatter(Xi_ref_1(1,:),Xi_ref_1(2,:),20,[1 0 0],'filled'); hold on
[hd] = scatter(Xi_ref_2(1,:),Xi_ref_2(2,:),20,[1 0 1],'filled'); hold on
[hs] = plot_ds_model(fig1, ds_lpv_2, att_g, limits,'medium'); hold on;
limits_ = limits + [-0.015 0.015 -0.015 0.015];
axis(limits_)
box on
grid on

xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

% plotGMMParameters( [Xi_ref_1 Xi_ref_2], est_labels_merged, merged_ds_gmm.Mu, merged_ds_gmm.Sigma, fig1);
title('Updated batch $b+1$ $\tilde{\theta}_{\gamma}^{b+1}$ and $\mathbf{\tilde{f}}(\xi)^{b+1}$', 'Interpreter','LaTex','FontSize',20)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Plot Choosen Lyapunov Function and derivative  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type of plot
contour = 1; % 0: surf, 1: contour
clear lyap_fun_comb lyap_der

P = eye(2);
% Lyapunov function
lyap_fun = @(x)lyapunov_function_PQLF(x, att_g, P);

% Derivative of Lyapunov function (gradV*f(x))
lyap_der = @(x)lyapunov_derivative_PQLF(x, att_g, P, ds_lpv_merged);
title_string_der = {'Lyapunov Function Derivative $\dot{V}(\xi)$'};

% if exist('h_lyap','var');     delete(h_lyap);     end
% if exist('h_lyap_der','var'); delete(h_lyap_der); end
h_lyap     = plot_lyap_fct(lyap_fun, contour, limits,  title_string, 0);
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);

h_lyap_der = plot_lyap_fct(lyap_der, contour, limits_,  title_string_der, 1);
[hd] = scatter(Xi_ref_1(1,:),Xi_ref_1(2,:),10,[1 1 0],'filled'); hold on
[hd] = scatter(Xi_ref_2(1,:),Xi_ref_2(2,:),10,[1 1 0],'filled'); hold on
xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
