function [lpv_rmse_train, lpv_edot_train, lpv_dtwd_train, lpv_rmse_test, lpv_edot_test, lpv_dtwd_test, lpv_K_s, names ] = extract_LPV_stats(lpv_stats, models_to_eval, opt)

lpv_rmse_train = []; lpv_edot_train = []; lpv_dtwd_train = [];
lpv_rmse_test = []; lpv_edot_test = [];   lpv_dtwd_test = []; 
names = [];ids = []; lpv_K_s = [];
for s=1:length(models_to_eval)
    
    s_ = models_to_eval(s);
    % Gathering Training stats
    lpv_rmse_train = [lpv_rmse_train lpv_stats{s_}.rmse_train(opt,:)'];    
    lpv_edot_train = [lpv_edot_train lpv_stats{s_}.edot_train(opt,:)'];
    dtwd_train = lpv_stats{s_}.dtwd_train(:,:,opt);
    lpv_dtwd_train = [lpv_dtwd_train dtwd_train(:)];

    
    % Gathering Testing stats
    lpv_rmse_test = [lpv_rmse_test lpv_stats{s_}.rmse_test(opt,:)'];
    lpv_edot_test = [lpv_edot_test lpv_stats{s_}.edot_test(opt,:)'];
    dtwd_test = lpv_stats{s_}.dtwd_test(:,:,opt);
    lpv_dtwd_test = [lpv_dtwd_test dtwd_test(:)];
    
    names{s} = lpv_stats{s_}.name;
    ids = [ids s_];
    lpv_K_s = [lpv_K_s lpv_stats{s_}.K_s'];
end

names

end