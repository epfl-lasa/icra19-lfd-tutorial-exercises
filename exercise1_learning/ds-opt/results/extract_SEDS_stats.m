function [seds_rmse_train, seds_edot_train, seds_dtwd_train, seds_rmse_test, seds_edot_test, seds_dtwd_test, seds_K_s, names ] = extract_SEDS_stats(seds_stats, models_to_eval)

seds_rmse_train = []; seds_edot_train = []; seds_dtwd_train = [];
seds_rmse_test = []; seds_edot_test = [];   seds_dtwd_test = []; 
names = [];ids = []; seds_K_s = [];
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
    
    names{s} = seds_stats{s_}.name;
    ids = [ids s_];
    seds_K_s = [seds_K_s seds_stats{s_}.K_s'];
end

names

end