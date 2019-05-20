%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluation of SEDS Models on LASA Dataset   %%
%  Train and compare SEDS models               %%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose DS LASA Dataset to load
clear all; close all; clc

% Select the motion models that you would like to evaluate
motion_model = [1:26];
seds_stats = [];

tic;
for m = 1:length(motion_model)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Load Chosen Dataset %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    show_trajectories = 0;
    sample = 2;
    [demos, name] = load_LASA_dataset_shape(motion_model(m));
    pause(0);
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Run SEDS Evaluation %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    repetitions = 10;
    rmse_train  = zeros(1,repetitions);
    rmse_test   = zeros(1,repetitions);
    edot_train  = zeros(1,repetitions);
    edot_test   = zeros(1,repetitions);
    dtwd_train  = zeros(4,repetitions);
    dtwd_test   = zeros(3,repetitions);
    K_s         = zeros(1,repetitions);
    
    % GMM Fitting Options
    est_options = [];
    est_options.type        = 1;   % GMM Estimation Alorithm Type
    est_options.maxK        = 15;  % Maximum Gaussians for Type 1/2
    est_options.do_plots    = 0;   % Plot Estimation Statistics
    est_options.fixed_K     = [];  % Fix K and estimate with EM
    est_options.exp_scaling = [];
    sample = 1;
    
    % SEDS Solver Options
    seds_options.tol_mat_bias = 10^-6; % A very small positive scalar to avoid
    seds_options.display = 1;          % An option to control whether the algorithm
    seds_options.tol_stopping=10^-9;   % A small positive scalar defining the stoppping
    seds_options.max_iter = 500;      % Maximum number of iteration for the solver [deseds_stats{m}fault: i_max=1000]
    seds_options.objective = 'mse';    % 'mse' or 'likelihood'
    
    for j = 1:repetitions
        
        % Do model selection to find adequate number of Gaussians
        [Priors0, Mu0, Sigma0] = discover_local_models([Xi_ref_train(:,1:sample:end); Xi_dot_ref_train(:,1:sample:end)], Xi_dot_ref_train(:,1:sample:end), est_options);
        nb_gaussians = length(Priors0);
        
        % Finding an initial guess for GMM's parameter
        [Priors_0, Mu_0, Sigma_0] = initialize_SEDS([Xi_ref_train(:,1:sample:end); Xi_dot_ref_train(:,1:sample:end)],nb_gaussians);
        
        % Run SEDS Solver
        [Priors, Mu, Sigma]= SEDS_Solver(Priors_0,Mu_0,Sigma_0,[Xi_ref_train(:,1:sample:end); Xi_dot_ref_train(:,1:sample:end)], seds_options);
        
        % Generate SEDS function
        clear ds_seds
        ds_seds = @(x) GMR_SEDS(Priors,Mu,Sigma,x,1:2,3:4);
        
        fprintf ('Repetition %d \n',j);
        % Compute Prediction RMSE on training data
        rmse_train(1,j) = mean(rmse_error(ds_seds, Xi_ref_train, Xi_dot_ref_train));
        fprintf('For %s dataset, got prediction RMSE on training set: %d \n', name, rmse_train(1,j));
        
        % Compute Prediction RMSE on testing data
        rmse_test(1,j) = mean(rmse_error(ds_seds, Xi_ref_test, Xi_dot_ref_test));
        fprintf('For %s dataset, got prediction RMSE on testing set: %d \n', name, rmse_test(1,j));
        
        % Compute Prediction e_dot on training data
        edot_train(1,j) = mean(edot_error(ds_seds, Xi_ref_train, Xi_dot_ref_train));
        fprintf('For %s dataset, got e_dot on training set: %d \n', name, edot_train(1,j));
        
        % Compute Prediction e_dot on testing data
        edot_test(1,j) = mean(edot_error(ds_seds, Xi_ref_test, Xi_dot_ref_test));
        fprintf('For %s dataset, got e_dot on testing set: %d \n', name, edot_test(1,j));
        
        % Record Parameters
        K_s (1,j) = nb_gaussians;
                
        % Simulating train/test trajectories to compute DTWD
        opt_sim = [];
        opt_sim.dt = 0.01;
        opt_sim.i_max = 3000;
        opt_sim.tol = 0.1;
        opt_sim.plot = 0;
        [x_seds_train, ~] = Simulation(x0_all_train ,[],ds_seds, opt_sim);
        [x_seds_test,  ~] = Simulation(x0_all_test ,[], ds_seds, opt_sim);
        
        % Compute Reproduction DTWD error on training data
        nb_traj_train       = size(x_seds_train,3);
        train_ref_traj_leng = size(Xi_ref_train,2)/nb_traj_train;
        for n=1:nb_traj_train
            dtwd_train(n,j) = dtw(x_seds_train(:,:,n)',Xi_ref_train(:,1+(n-1)*train_ref_traj_leng:n*train_ref_traj_leng)')';
        end
        fprintf('For %s dataset, got reproduction DTWD on training set: %2.4f +/- %2.4f \n', name, mean(dtwd_train(:,j)),std(dtwd_train(:,j)));
    
        % Compute Reproduction DTWD error on testing data
        nb_traj_test       = size(x_seds_test,3);
        test_ref_traj_leng = size(Xi_ref_test,2)/nb_traj_test;
        for n=1:nb_traj_test
            dtwd_test(n,j) = dtw(x_seds_test(:,:,n)',Xi_ref_test(:,1+(n-1)*test_ref_traj_leng:n*test_ref_traj_leng)')';
        end
        fprintf('For %s dataset, got reproduction DTWD on testing set: %2.4f +/- %2.4f \n', name, mean(dtwd_test(:,j)),std(dtwd_test(:,j)));
    
    end
    
    seds_stats{m}.motion_model = motion_model(m);
    seds_stats{m}.name = name;
    seds_stats{m}.K_s = K_s;
    seds_stats{m}.rmse_train = rmse_train;
    seds_stats{m}.rmse_test  = rmse_test;
    seds_stats{m}.edot_train = edot_train;
    seds_stats{m}.edot_test  = edot_test;
    seds_stats{m}.dtwd_train = dtwd_train;
    seds_stats{m}.dtwd_test  = dtwd_test;
    
end
toc;

%% Compute and Plot Results
seds_rmse_train = []; seds_edot_train = []; seds_dtwd_train = [];
seds_rmse_test = []; seds_edot_test = [];   seds_dtwd_test = []; 
names = [];ids = []; seds_K_s = [];

% Options
models_to_eval = [6];
show_plot = 0;
for s=1:length(models_to_eval)
    
    s_ = models_to_eval(s)';
    % Gathering Training stats
    seds_rsme_model_ = seds_stats{s_}.rmse_train';
    seds_rmse_train = [seds_rmse_train seds_rsme_model_];
    
    seds_edot_model = seds_stats{s_}.edot_train';
    seds_edot_train = [seds_edot_train seds_edot_model];
    
    seds_dtwd_model = seds_stats{s_}.dtwd_train;
    seds_dtwd_train = [seds_dtwd_train seds_dtwd_model(:)];
    
    % Gathering Testing stats
    seds_rsme_model = [seds_stats{s_}.rmse_test]';
    seds_rmse_test  = [seds_rmse_test seds_rsme_model];
    
    seds_edot_model = [seds_stats{s_}.edot_test]';
    seds_edot_test = [seds_edot_test seds_edot_model];
    
    seds_dtwd_model = seds_stats{s_}.dtwd_test;
    seds_dtwd_test = [seds_dtwd_test seds_dtwd_model(:)];
    
    names{s} = seds_stats{s_}.name
    ids = [ids s_];
    seds_K_s = [seds_K_s seds_stats{s_}.K_s'];
end
if show_plot
    figure('Color',[1 1 1]);
    subplot(3,2,1)
    % notBoxPlot(seds_rmse_test,ids,'style','patch')
    violinplot(seds_rmse_train,ids);
    grid on; box on;
    title('SEDS Performance on Training Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('Motion Models','Interpreter','LaTex', 'FontSize',20);
    
    subplot(3,2,3)
    % notBoxPlot(seds_edot_test,ids,'style','patch')
    violinplot(seds_edot_train,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    subplot(3,2,2)
    % notBoxPlot(seds_rmse_test,ids,'style','patch')
    violinplot(seds_rmse_test,ids);
    grid on; box on;
    title('SEDS Performance on Testing Sets','Interpreter','LaTex', 'FontSize',20);
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20);
    xlabel('Motion Models','Interpreter','LaTex', 'FontSize',20);
    
    subplot(3,2,4)
    % notBoxPlot(seds_edot_test,ids,'style','patch')
    violinplot(seds_edot_test,ids);
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20);
    
    
    subplot(3,2,5)
    % notBoxPlot(seds_edot_test,ids,'style','patch')
    violinplot(seds_dtwd_train,ids);
    grid on; box on;
    ylabel('DTWD','Interpreter','LaTex', 'FontSize',20);
    
    subplot(3,2,6)
    % notBoxPlot(seds_edot_test,ids,'style','patch')
    violinplot(seds_dtwd_test,ids);
    grid on; box on;
    ylabel('DTWD','Interpreter','LaTex', 'FontSize',20);
    
end
figure('Color',[1 1 1]);
notBoxPlot(seds_K_s,ids,'style','sdline');
% boxplot(K_s,ids);
grid on; box on;
title('Optimal $K$ Choosen via Model Selection for SEDS $p(\xi,\dot{\xi})$','Interpreter','LaTex', 'FontSize',20);
ylabel('K','Interpreter','LaTex', 'FontSize',20);
xlabel('Motion Models','Interpreter','LaTex', 'FontSize',20);
fprintf('Optimal K %2.2f +/- %2.4f\n', mean(seds_K_s), std(seds_K_s));

%% Overall Performance
% clc;
fprintf('**************** Performance of SEDS ****************\n');
seds_train_mean = mean(seds_rmse_train(:));
seds_train_std  = std(seds_rmse_train(:));
fprintf('Overall SEDS RMSE on Training Set %2.4f +/- %2.4f\n', seds_train_mean, seds_train_std);
seds_test_mean = mean(seds_rmse_test(:));
seds_test_std  = std(seds_rmse_test(:));
fprintf('Overall SEDS RMSE on Testing Set %2.4f +/- %2.4f\n', seds_test_mean, seds_test_std);

seds_train_mean_edot = mean(seds_edot_train(:));
seds_train_std_edot  = std(seds_edot_train(:));
fprintf('Overall SEDS e-dot on Training Set %2.4f +/- %2.4f\n', seds_train_mean_edot, seds_train_std_edot);
seds_test_mean_edot = mean(seds_edot_test(:));
seds_test_std_edot  = std(seds_edot_test(:));
fprintf('Overall SEDS e-dot on Testing Set %2.4f +/- %2.4f\n', seds_test_mean_edot, seds_test_std_edot);


seds_train_mean_dtwd = mean(seds_dtwd_train(:));
seds_train_std_dtwd  = std(seds_dtwd_train(:));
fprintf('Overall SEDS DTWD on Training Set %2.4f +/- %2.4f\n', seds_train_mean_dtwd, seds_train_std_dtwd);
seds_test_mean_dtwd = mean(seds_dtwd_test(:));
seds_test_std_dtwd  = std(seds_dtwd_test(:));
fprintf('Overall SEDS DTWD on Testing Set %2.4f +/- %2.4f\n', seds_test_mean_dtwd, seds_test_std_dtwd);


figure('Color',[1 1 1]);
subplot(3,1,1);
grid on; box on;
violinplot([seds_rmse_train(:) seds_rmse_test(:)],{'SEDS-Train','SEDS-Test'});
ylabel('RMSE','Interpreter','LaTex', 'FontSize',20)
title('Overall Reproduction Accuracy', 'Interpreter','LaTex', 'FontSize',20)

subplot(3,1,2);
violinplot([seds_edot_train(:) seds_edot_test(:)],{'SEDS-Train','SEDS-Test'});
grid on; box on;
ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20)

subplot(3,1,3);
dummy_test = seds_dtwd_test(:);
dummy_test = [dummy_test; mean(seds_dtwd_test(:))*ones(length(seds_dtwd_train(:)) - length(seds_dtwd_test(:)),1)];
violinplot([seds_dtwd_train(:)  dummy_test],{'SEDS-Train','SEDS-Test'}); hold on;

grid on; box on;
ylabel('DTWD','Interpreter','LaTex', 'FontSize',20)