%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Comparison Script on LASA Motion Library Dataset      %%
%       Load stats from different algorithms and compare!      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load SEDS stats
clear all; close all; clc;
load('seds_stats_all');

%% Load EM stats
em_lpv_stats = [];
for i=1:26
   tmp          = load(sprintf('./em-per-motion/em_lpv_stats_%d.mat', i));
   struct_name  = sprintf('lpv_stats_%d',i);
   em_lpv_stats{i} = tmp.(struct_name){1};
end

%% Load PC stats
pc_lpv_stats = [];
for i=1:26
   tmp          = load(sprintf('./pc-per-motion/pc_lpv_stats_%d.mat', i));
   struct_name  = sprintf('lpv_stats_%d',i);
   pc_lpv_stats{i} = tmp.(struct_name){1};
end

%% Models to Evaluate
% For individual motions
models_to_eval = [6]; 

% Challenging motions for SEDS
% models_to_eval = [2,3,7,8,10,11,14,18,20,22,26]; 

% For full LASA Dataset (excluing multi-models)
% models_to_eval = [1:26]; 

% Extract stats from SEDS Solver
[seds_rmse_train, seds_edot_train, seds_dtwd_train, ...
 seds_rmse_test, seds_edot_test, seds_dtwd_test, seds_K_s, names ] = extract_SEDS_stats(seds_stats, models_to_eval);

% Extract stats from EM-LPV O1-3 Solver
em_lpv_rmse_train = []; em_lpv_edot_train = []; em_lpv_dtwd_train = [];
em_lpv_rmse_test  = []; em_lpv_edot_test = []; em_lpv_dtwd_test = []; em_lpv_K_s = [];
for opt=1:3
    [em_lpv_rmse_train(:,:,opt), em_lpv_edot_train(:,:,opt), em_lpv_dtwd_train(:,:,opt), ...
     em_lpv_rmse_test(:,:,opt),  em_lpv_edot_test(:,:,opt), em_lpv_dtwd_test(:,:,opt), em_lpv_K_s(:,:,opt), ~ ] = extract_LPV_stats(em_lpv_stats, models_to_eval, opt);
end

% Extract stats from EM-LPV O1-3 Solver
pc_lpv_rmse_train = []; pc_lpv_edot_train = []; pc_lpv_dtwd_train = [];
pc_lpv_rmse_test  = []; pc_lpv_edot_test = [];  pc_lpv_dtwd_test = []; pc_lpv_K_s = [];
for opt=1:3
    [pc_lpv_rmse_train(:,:,opt), pc_lpv_edot_train(:,:,opt), pc_lpv_dtwd_train(:,:,opt), ...
     pc_lpv_rmse_test(:,:,opt), pc_lpv_edot_test(:,:,opt), pc_lpv_dtwd_test(:,:,opt), pc_lpv_K_s(:,:,opt), ~ ] = extract_LPV_stats(pc_lpv_stats, models_to_eval, opt);
end

%% Plot individual motion performance
clc;
computeMetricsPerMotion(seds_rmse_train, seds_edot_train, seds_dtwd_train, ...
                       em_lpv_rmse_train, em_lpv_edot_train, em_lpv_dtwd_train, ...
                       pc_lpv_rmse_train, pc_lpv_edot_train, pc_lpv_dtwd_train, names, 'train', 1);
                   

computeMetricsPerMotion(seds_rmse_test, seds_edot_test, seds_dtwd_test, ...
                       em_lpv_rmse_test, em_lpv_edot_test, em_lpv_dtwd_test, ...
                       pc_lpv_rmse_test, pc_lpv_edot_test, pc_lpv_dtwd_test, names, 'test',1);                   


%% Plot overall motion performance
clc;
computeMetricsOverall(seds_rmse_train, seds_edot_train, seds_dtwd_train, ...
                      em_lpv_rmse_train, em_lpv_edot_train, em_lpv_dtwd_train, ...
                      pc_lpv_rmse_train, pc_lpv_edot_train, pc_lpv_dtwd_train, 'train',1);

computeMetricsOverall(seds_rmse_test, seds_edot_test, seds_dtwd_test, ...
                      em_lpv_rmse_test, em_lpv_edot_test, em_lpv_dtwd_test, ...
                      pc_lpv_rmse_test, pc_lpv_edot_test, pc_lpv_dtwd_test, 'test',1);
