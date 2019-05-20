%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluation of Learning Schemes for LPVDS Models on LASA Dataset    %%
%  Train and compare a series of LPVDS models with different GMM      %%   
%  fitting approaches and optimization variants                       %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      LOAD Motion from LASA DATASET (30 motions)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
clear all; close all; clc

% Select the motion models that you would like to evaluate
motion_model = [1];
seds_stats = [];

tic;
for m = 1:length(motion_model)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Load Chosen Dataset %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    sample = 2;
    show_trajectories = 0;
    [demos, name] = load_LASA_dataset_shape(m);
    
    % Global Attractor of DS
    att_g = [0 0]';
    
    % Generate Training Data
    Data_train = []; x0_all_train = [];
    for l=1:4
        % Check where demos end and shift
        data_ = [demos{l}.pos(:,1:sample:end); demos{l}.vel(:,1:sample:end);];
        Data_train = [Data_train data_];
        x0_all_train = [x0_all_train data_(1:2,20)];
        clear data_
    end
    % Position/Velocity Trajectories
    Xi_ref_train     = Data_train(1:2,:);
    Xi_dot_ref_train = Data_train(3:end,:);
    
    % Generate Testing Data
    Data_test = []; x0_all_test = [];
    for l=5:7
        % Check where demos end and shift
        data_ = [demos{l}.pos(:,1:sample:end); demos{l}.vel(:,1:sample:end);];
        Data_test = [Data_test data_];
        x0_all_test = [x0_all_test data_(1:2,20)];
        clear data_
    end
    % Position/Velocity Trajectories
    Xi_ref_test     = Data_test(1:2,:);
    Xi_dot_ref_test = Data_test(3:end,:);
    
    if show_trajectories
        figure('Color',[1 1 1]);
        scatter(Xi_ref_train(1,:),Xi_ref_train(2,:),10,[1 0 0],'filled'); hold on
        scatter(Xi_ref_test(1,:), Xi_ref_test(2,:),10,[0 0 1],'filled'); hold on
        grid on; box on;
        xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
        ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
        legend('Training Data', 'Testing Data')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Run LPVDS Evaluation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%% GMM Estimation Algorithm %%%%%%%%%%%%%%%%%%%
    % 0: Physically-Consistent Non-Parametric (Collapsed Gibbs Sampler)
    % 1: GMM-EM Model Selection via BIC
    % 2: GMM via Competitive-EM
    % 3: CRP-GMM via Collapsed Gibbs Sampler
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gmm_fit = 1;
    
    %%%%%%%%%%%%%%%%%%% DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%
    % Type of constraints/optimization
    % constr_type = 0;      % 0:'convex':     A' + A < 0
    % 1:'non-convex': A'P + PA < 0
    % 2:'non-convex': A'P + PA < -Q given P
    % init_cvx    = 0;      % 0/1: initialize non-cvx problem with cvx
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constr_type = [0,1,2];
    
    % Repetitions
    repetitions = 10;
    
    % Errors and parameters
    rmse_train  = zeros(length(constr_type),repetitions);
    rmse_test   = zeros(length(constr_type),repetitions);
    edot_train  = zeros(length(constr_type),repetitions);
    edot_test   = zeros(length(constr_type),repetitions);
    K_s         = zeros(1,repetitions);
    
    
    % GMM Fitting Options
    est_options = [];
    est_options.type        = gmm_fit;    % GMM Estimation Alorithm Type
    est_options.maxK        = 15;         % Maximum Gaussians for Type 1/2
    est_options.do_plots    = 0;          % Plot Estimation Statistics
    est_options.adjusts_C   = 0;          % Adjust Sigmas
    est_options.fixed_K     = [];         % Fix K and estimate with EM
    est_options.exp_scaling = 1;          % Scaling for the similarity to improve locality
    sample = 2;
    
    for j = 1:repetitions
        
        % Discover Local Models
        [Priors0, Mu0, Sigma0] = discover_local_models(Xi_ref_train(:,1:sample:end), Xi_ref_train(:,1:sample:end), est_options);
        
        % Create GMM Structure
        clear ds_gmm; ds_gmm.Mu = Mu0; ds_gmm.Sigma = Sigma0;
        ds_gmm.Priors = Priors0;
        est_K      = length(Priors0);
        
        if est_options.adjusts_C  == 1
            tot_scale_fact = 1; rel_scale_fact = 0.15;
            Sigma = adjust_Covariances(Sigma0, tot_scale_fact, rel_scale_fact);
            ds_gmm.Sigma;
        end
        
        % Record Parameters
        K_s (1,j) = est_K;
        
        for i = 1:length(constr_type)
            %%%%%%%%%%%%%%%%%%% DS OPTIMIZATION OPTIONS %%%%%%%%%%%%%%%%%%%%
            % Type of constraints/optimization
            % constr_type = 0;      % 0:'convex':     A' + A < 0
            % 1:'non-convex': A'P + PA < 0
            % 2:'non-convex': A'P + PA < -Q given P
            % init_cvx    = 0;      % 0/1: initialize non-cvx problem with cvx
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%  LPV system sum_{k=1}^{K}\gamma_k(xi)(A_kxi + b_k) %%%%%%%%
            if constr_type(i) == 0 && constr_type(i) == 1 %% This should be ||
                init_cvx    = 0; P_opt = [];
            else
                [Vxf] = learn_wsaqf(Data_train, att_g);
                P_opt = Vxf.P(:,:,1);
                init_cvx    = 0;
            end
            
            [A_g, b_g, P_est] = optimize_lpv_ds_from_data(Data_train, att_g, constr_type(i), ds_gmm, P_opt, init_cvx);
            
            % Create DS function handle
            clear ds_lpv
            ds_lpv = @(x) lpv_ds(x, ds_gmm, A_g, b_g);
            
            fprintf ('Repetition %d \n',j);
            % Compute RMSE on training data
            rmse_train(i,j) = mean(rmse_error(ds_lpv, Xi_ref_train, Xi_dot_ref_train));
            fprintf('For %s dataset with (O%d), got prediction RMSE on training set: %d \n', name, constr_type(i)+1, rmse_train(i,j));
            
            % Compute RMSE on testing data
            rmse_test(i,j) = mean(rmse_error(ds_lpv, Xi_ref_test, Xi_dot_ref_test));
            fprintf('For %s dataset with (O%d), got prediction RMSE on testing set: %d \n', name, constr_type(i)+1, rmse_test(i,j));
            
            % Compute e_dot on training data
            edot_train(i,j) = mean(edot_error(ds_lpv, Xi_ref_train, Xi_dot_ref_train));
            fprintf('For %s dataset with (O%d), got e_dot on training set: %d \n', name, constr_type(i)+1, edot_train(i,j));
            
            % Compute e_dot on testing data
            edot_test(i,j) = mean(edot_error(ds_lpv, Xi_ref_test, Xi_dot_ref_test));
            fprintf('For %s dataset with (O%d), got e_dot on testing set: %d \n', name, constr_type(i)+1, edot_test(i,j));

        end
    end
    lpv_stats{m}.motion_model = motion_model(m);
    lpv_stats{m}.name = name;
    lpv_stats{m}.K_s = K_s;
    lpv_stats{m}.rmse_train = rmse_train;
    lpv_stats{m}.rmse_test  = rmse_test;
    lpv_stats{m}.edot_train = edot_train;
    lpv_stats{m}.edot_test  = edot_test;
end
toc;

%% Compute and Plot Results
show_plot = 1;

lpv_rmse_train_O1 = []; lpv_edot_train_O1 = [];
lpv_rmse_test_O1 = []; lpv_edot_test_O1 = [];

lpv_rmse_train_O2 = []; lpv_edot_train_O2 = [];
lpv_rmse_test_O2 = []; lpv_edot_test_O2 = [];

lpv_rmse_train_O3 = []; lpv_edot_train_O3 = [];
lpv_rmse_test_O3 = []; lpv_edot_test_O3 = [];

names = [];ids = []; K_s = [];
% num_models = length(lpv_stats);
num_models = 26;
for s=1:num_models
    % Gathering Training stats
    lpv_rmse_train_O1 = [lpv_rmse_train_O1 lpv_stats{s}.rmse_train(1,:)'];
    lpv_rmse_train_O2 = [lpv_rmse_train_O2 lpv_stats{s}.rmse_train(2,:)'];
    lpv_rmse_train_O3 = [lpv_rmse_train_O3 lpv_stats{s}.rmse_train(3,:)'];
    
    lpv_edot_train_O1 = [lpv_edot_train_O1 lpv_stats{s}.edot_train(1,:)'];
    lpv_edot_train_O2 = [lpv_edot_train_O2 lpv_stats{s}.edot_train(2,:)'];
    lpv_edot_train_O3 = [lpv_edot_train_O3 lpv_stats{s}.edot_train(3,:)'];
    
    % Gathering Testing stats
    lpv_rmse_test_O1 = [lpv_rmse_test_O1 lpv_stats{s}.rmse_test(1,:)'];
    lpv_rmse_test_O2 = [lpv_rmse_test_O2 lpv_stats{s}.rmse_test(2,:)'];
    lpv_rmse_test_O3 = [lpv_rmse_test_O3 lpv_stats{s}.rmse_test(3,:)'];
    
    lpv_edot_test_O1 = [lpv_edot_test_O1 lpv_stats{s}.edot_test(1,:)'];
    lpv_edot_test_O2 = [lpv_edot_test_O2 lpv_stats{s}.edot_test(2,:)'];
    lpv_edot_test_O3 = [lpv_edot_test_O3 lpv_stats{s}.edot_test(3,:)'];
    
    names{s} = lpv_stats{s}.name;
    ids = [ids s];
    K_s = [K_s lpv_stats{s}.K_s'];
end
if show_plot
    %%%%% Plot Performance of (O1)
    figure('Color',[1 1 1]);
    subplot(2,2,1)
    violinplot(lpv_rmse_train_O1,ids)
%     set(gca,'xtick',[1:num_models],'xticklabel', names)
    % xtickangle(45)
    grid on; box on;
    title('PC-GMM LPV-DS (O1) Performance on Training Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,3)
    violinplot(lpv_edot_train_O1,ids)
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,2)
    violinplot(lpv_rmse_test_O1,ids)
    grid on; box on;
    title('PC-GMM LPV-DS (O1) Performance on Testing Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,4)
    violinplot(lpv_edot_test_O1,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    %%%%% Plot Performance of (O2)
    figure('Color',[1 1 1]);
    subplot(2,2,1)
    violinplot(lpv_rmse_train_O2,ids);
    grid on; box on;
    title('PC-GMM LPV-DS (O2) Performance on Training Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,3)
    violinplot(lpv_edot_train_O2,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,2)
    violinplot(lpv_rmse_test_O2,ids);
    grid on; box on;
    title('PC-GMM LPV-DS (O2) Performance on Testing Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,4)
    violinplot(lpv_edot_test_O2,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    
    %%%%% Plot Performance of (O3)
    figure('Color',[1 1 1]);
    subplot(2,2,1)
    violinplot(lpv_rmse_train_O3,ids);
    grid on; box on;
    title('PC-GMM LPV-DS (O3) Performance on Training Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,3)
    violinplot(lpv_edot_train_O3,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,2)
    violinplot(lpv_rmse_test_O3,ids);
    grid on; box on;
    title('PC-GMM LPV-DS (O3) Performance on Testing Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);
    
    subplot(2,2,4)
    violinplot(lpv_edot_test_O3,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);    
end

% K's
figure('Color',[1 1 1]);
notBoxPlot(K_s,ids,'style','patch');
grid on; box on;
title('Optimal $K$ via PC-GMM','Interpreter','LaTex', 'FontSize',20);
ylabel('K','Interpreter','LaTex', 'FontSize',20);
xlabel('LASA Motion Model Library','Interpreter','LaTex', 'FontSize',20);


%% Overall Performance
% clc;
show_plot = 1;
fprintf('\n**************** Performance of PC-GMM LPV-DS (O1) ****************\n');
lpv_train_mean_O1 = mean(lpv_rmse_train_O1(:));
lpv_train_std_O1  = std(lpv_rmse_train_O1(:));
fprintf('Overall RMSE on Training Set %2.2f +/- %2.2f\n', lpv_train_mean_O1, lpv_train_std_O1);
lpv_test_mean_O1 = mean(lpv_rmse_test_O1(:));
lpv_test_std_O1  = std(lpv_rmse_test_O1(:));
fprintf('Overall RMSE on Testing Set %2.2f +/- %2.2f\n', lpv_test_mean_O1, lpv_test_std_O1);

lpv_train_mean_edot_O1 = mean(lpv_edot_train_O1(:));
lpv_train_std_edot_O1  = std(lpv_edot_train_O1(:));
fprintf('Overall e-dot on Training Set %2.2f +/- %2.2f\n', lpv_train_mean_edot_O1, lpv_train_std_edot_O1);
lpv_test_mean_edot_O1 = mean(lpv_edot_test_O1(:));
lpv_test_std_edot_O1  = std(lpv_edot_test_O1(:));
fprintf('Overall e-dot on Testing Set %2.2f +/- %2.2f\n', lpv_test_mean_edot_O1, lpv_test_std_edot_O1);

fprintf('**************** Performance of PC-GMM LPV-DS (O2) ****************\n');
lpv_train_mean_O2 = mean(lpv_rmse_train_O2(:));
lpv_train_std_O2  = std(lpv_rmse_train_O2(:));
fprintf('Overall RMSE on Training Set %2.2f +/- %2.2f\n', lpv_train_mean_O2, lpv_train_std_O2);
lpv_test_mean_O2 = mean(lpv_rmse_test_O2(:));
lpv_test_std_O2  = std(lpv_rmse_test_O2(:));
fprintf('Overall RMSE on Testing Set %2.2f +/- %2.2f\n', lpv_test_meanDave Chen Pao_O2, lpv_test_std_O2);

lpv_train_mean_edot_O2 = mean(lpv_edot_train_O2(:));
lpv_train_std_edot_O2  = std(lpv_edot_train_O2(:));
fprintf('Overall e-dot on Training Set %2.2f +/- %2.2f\n', lpv_train_mean_edot_O2, lpv_train_std_edot_O2);
lpv_test_mean_edot_O2 = mean(lpv_edot_test_O2(:));
lpv_test_std_edot_O2  = std(lpv_edot_test_O2(:));
fprintf('Overall e-dot on Testing Set %2.2f +/- %2.2f\n', lpv_test_mean_edot_O2, lpv_test_std_edot_O2);

fprintf('**************** Performance of PC-GMM LPV-DS (O3) ****************\n');
lpv_train_mean_O3 = mean(lpv_rmse_train_O3(:));
lpv_train_std_O3  = std(lpv_rmse_train_O3(:));
fprintf('Overall RMSE on Training Set %2.2f +/- %2.2f\n', lpv_train_mean_O3, lpv_train_std_O3);
lpv_test_mean_O3 = mean(lpv_rmse_test_O3(:));
lpv_test_std_O3  = std(lpv_rmse_test_O3(:));
fprintf('Overall RMSE on Testing Set %2.2f +/- %2.2f\n', lpv_test_mean_O3, lpv_test_std_O3);

lpv_train_mean_edot_O3 = mean(lpv_edot_train_O3(:));
lpv_train_std_edot_O3  = std(lpv_edot_train_O3(:));
fprintf('Overall e-dot on Training Set %2.2f +/- %2.2f\n', lpv_train_mean_edot_O3, lpv_train_std_edot_O3);
lpv_test_mean_edot_O3 = mean(lpv_edot_test_O3(:));
lpv_test_std_edot_O3  = std(lpv_edot_test_O3(:));
fprintf('Overall e-dot on Testing Set %2.2f +/- %2.2f\n', lpv_test_mean_edot_O3, lpv_test_std_edot_O3);

if show_plot
    figure('Color',[1 1 1]);
    subplot(2,2,1);
    grid on; box on;
    violinplot([lpv_rmse_train_O1(:) lpv_rmse_train_O2(:) lpv_rmse_train_O3(:)],{'(O1)','(O2)','(O3)'});
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20)
    title('Overall Reproduction Accuracy PC-GMM LPV-DS (Training) ', 'Interpreter','LaTex', 'FontSize',18)
    
    subplot(2,2,3);
    violinplot([lpv_edot_train_O1(:) lpv_edot_train_O2(:) lpv_edot_train_O3(:)],{'(O1)','(O2)','(O3)'});
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20)
    
    subplot(2,2,2);
    grid on; box on;
    violinplot([lpv_rmse_test_O1(:) lpv_rmse_test_O2(:) lpv_rmse_test_O3(:)],{'(O1)','(O2)','(O3)'});
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20)
    title('Overall Reproduction Accuracy PC-GMM LPV-DS (Testing) ', 'Interpreter','LaTex', 'FontSize',18)
    
    subplot(2,2,4);
    violinplot([lpv_edot_test_O1(:) lpv_edot_test_O2(:) lpv_edot_test_O3(:)],{'(O1)','(O2)','(O3)'});
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20)
        
end

