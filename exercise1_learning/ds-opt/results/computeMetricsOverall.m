function [] =  computeMetricsOverall(seds_rmse, seds_edot, seds_dtwd, em_lpv_rmse, em_lpv_edot, em_lpv_dtwd, pc_lpv_rmse, pc_lpv_edot, pc_lpv_dtwd, data_type, face_color)

fprintf('**************** Performance of SEDS (Overall %s-set) ****************\n',data_type);
seds_rmse_mean = mean(seds_rmse(:));
seds_rmse_std  = std(seds_rmse(:));
fprintf('Overall SEDS RMSE on %s Set %2.4f +/- %2.4f\n', data_type, seds_rmse_mean, seds_rmse_std);

seds_edot_mean = mean(seds_edot(:));
seds_edot_std  = std(seds_edot(:));
fprintf('Overall SEDS e-dot on %s Set %2.4f +/- %2.4f\n', data_type, seds_edot_mean, seds_edot_std);

seds_dtwd_mean = mean(seds_dtwd(:));
seds_dtwd_std  = std(seds_dtwd(:));
fprintf('Overall SEDS DTWD on %s Set %2.4f +/- %2.4f\n', data_type, seds_dtwd_mean, seds_dtwd_std);

opt_em = size(em_lpv_rmse,3);
for opt=1:opt_em
    em_lpv_rmse_all = em_lpv_rmse(:,:,opt);
    em_lpv_rmse_mean(opt) = mean(em_lpv_rmse_all(:));
    em_lpv_rmse_std(opt)  = std(em_lpv_rmse_all(:));
    fprintf('EM-GMM LPV (O%d) RMSE on %s-set %2.4f +/- %2.4f\n', opt, data_type, em_lpv_rmse_mean(opt), em_lpv_rmse_std(opt));
        
    em_lpv_edot_all = em_lpv_edot(:,:,opt);
    em_lpv_edot_mean(opt) = mean(em_lpv_edot_all(:));
    em_lpv_edot_std(opt)  = std(em_lpv_edot_all(:));
    fprintf('EM-GMM LPV (O%d) e-dot on %s-set %2.4f +/- %2.4f\n', opt, data_type, em_lpv_edot_mean(opt), em_lpv_edot_std(opt));
    
    em_lpv_dtwd_all = em_lpv_dtwd(:,:,opt);
    em_lpv_dtwd_mean(opt) = mean(em_lpv_dtwd_all(:));
    em_lpv_dtwd_std(opt)  = std(em_lpv_dtwd_all(:));
    fprintf('EM-GMM LPV (O%d) DTWD on %s-set %2.4f +/- %2.4f\n', opt, data_type, em_lpv_dtwd_mean(opt), em_lpv_dtwd_std(opt));
end

opt_em = size(em_lpv_rmse,3);
for opt=1:opt_em
    pc_lpv_rmse_all = pc_lpv_rmse(:,:,opt);
    pc_lpv_rmse_mean(opt) = mean(pc_lpv_rmse_all(:));
    pc_lpv_rmse_std(opt)  = std(pc_lpv_rmse_all(:));
    fprintf('PC-GMM LPV (O%d) RMSE on %s-set %2.4f +/- %2.4f\n', opt, data_type, pc_lpv_rmse_mean(opt), pc_lpv_rmse_std(opt));
        
    pc_lpv_edot_all = pc_lpv_edot(:,:,opt);
    pc_lpv_edot_mean(opt) = mean(pc_lpv_edot_all(:));
    pc_lpv_edot_std(opt)  = std(pc_lpv_edot_all(:));
    fprintf('PC-GMM LPV (O%d) e-dot on %s-set %2.4f +/- %2.4f\n', opt, data_type, pc_lpv_edot_mean(opt), pc_lpv_edot_std(opt));
    
    pc_lpv_dtwd_all = pc_lpv_dtwd(:,:,opt);
    pc_lpv_dtwd_mean(opt) = mean(pc_lpv_dtwd_all(:));
    pc_lpv_dtwd_std(opt)  = std(pc_lpv_dtwd_all(:));
    fprintf('PC-GMM LPV (O%d) DTWD on %s-set %2.4f +/- %2.4f\n', opt, data_type, pc_lpv_dtwd_mean(opt), pc_lpv_dtwd_std(opt));
end


% Plot statistics
figure('Color',[1 1 1]);
subplot(3,1,1);
face_C = [0.75 0.75 1;
    1 0.75 0.75;
    1 0.5 0.5;
    1 0.25 0.25;
    0.75 1 0.75;
    0.5 1 0.5;
    0 1 0.25;
    ];
edge_C = [0  0    0.9;
    0.9 0 0;
    0.9 0 0;
    0.9 0 0;
    0   0.85 0;
    0   0.85 0;
    0   0.85 0;
    ];

% edge_C = [rgb('MidNightBlue');
%     rgb('Indigo');
%     rgb('DarkSlateBlue');
%     rgb('BlueViolet');
%     rgb('RoyalBlue');
%     rgb('DodgerBlue');
%     rgb('SteelBlue');
%     ];
% 
% face_C = [rgb('MediumBlue');
%     rgb('DarkMagenta');
%     rgb('DarkOrchid');
%     rgb('DarkViolet');
%     rgb('DeepSkyBlue');
%     rgb('SkyBlue');
%     rgb('LightSkyBlue');
%     ];

names_ = {'S','E(O1)','E(O2)','E(O3)', 'PC(O1)', 'PC(O2)', 'PC(O3)'};
X = [1 1.25 1.5 1.75 2 2.25 2.5];

title(sprintf('Overall ($%s$) Performance on LASA Library',data_type), 'Interpreter','LaTex', 'FontSize',16)
Y_rmse = [seds_rmse_mean em_lpv_rmse_mean(1) em_lpv_rmse_mean(2) em_lpv_rmse_mean(3) ...
          pc_lpv_rmse_mean(1) pc_lpv_rmse_mean(2) pc_lpv_rmse_mean(3)];
E_rmse = 2*[seds_rmse_std em_lpv_rmse_std(1) em_lpv_rmse_std(2) em_lpv_rmse_std(3) ...
          pc_lpv_rmse_std(1) pc_lpv_rmse_std(2) pc_lpv_rmse_std(3)];
if face_color
    h = superbar(X, Y_rmse, 'BarFaceColor', face_C,'BarEdgeColor', edge_C, 'E', E_rmse, 'ErrorbarStyle', 'T');
    set(h, 'LineWidth', 1);
else
    h = superbar(X, Y_rmse, 'BarFaceColor', 'none','BarEdgeColor', edge_C, 'E', E_rmse, 'ErrorbarStyle', 'T');
    set(h, 'LineWidth', 2);
end
xlim([0.85 2.65]);
grid on; box on;
ylabel('RMSE','Interpreter','LaTex', 'FontSize',16)
set(gca,'xtick',[1:0.25:2.5],'xticklabel',names_)

subplot(3,1,2);
Y_edot = [seds_edot_mean em_lpv_edot_mean(1) em_lpv_edot_mean(2) em_lpv_edot_mean(3) ...
          pc_lpv_edot_mean(1) pc_lpv_edot_mean(2) pc_lpv_edot_mean(3)];
E_edot = 2*[seds_edot_std em_lpv_edot_std(1) em_lpv_edot_std(2) em_lpv_edot_std(3) ...
          pc_lpv_edot_std(1) pc_lpv_edot_std(2) pc_lpv_edot_std(3)];
if face_color
    h = superbar(X, Y_edot, 'BarFaceColor', face_C,'BarEdgeColor', edge_C,'E', E_edot, 'ErrorbarStyle', 'T');
    set(h, 'LineWidth', 1);
else
    h = superbar(X, Y_edot, 'BarFaceColor', 'none','BarEdgeColor', edge_C, 'E', E_edot, 'ErrorbarStyle', 'T');
    set(h, 'LineWidth', 2);
end
xlim([0.85 2.65]);
grid on; box on;
ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',16)
set(gca,'xtick',[1:0.25:2.5],'xticklabel',names_)

subplot(3,1,3);
Y_dtwd = [seds_dtwd_mean em_lpv_dtwd_mean(1) em_lpv_dtwd_mean(2) em_lpv_dtwd_mean(3) ...
    pc_lpv_dtwd_mean(1) pc_lpv_dtwd_mean(2) pc_lpv_dtwd_mean(3)];
E_dtwd = 2*[seds_dtwd_std em_lpv_dtwd_std(1) em_lpv_dtwd_std(2) em_lpv_dtwd_std(3) ...
    pc_lpv_dtwd_std(1) pc_lpv_dtwd_std(2) pc_lpv_dtwd_std(3)];
subplot(3,1,3);
if face_color
    h = superbar(X, Y_dtwd, 'BarFaceColor',face_C,'BarEdgeColor', edge_C, 'E', E_dtwd, 'ErrorbarStyle', 'T');
    set(h, 'LineWidth', 1);
else
    h = superbar(X, Y_dtwd, 'BarFaceColor','none','BarEdgeColor', edge_C, 'E', E_dtwd, 'ErrorbarStyle', 'T');
    set(h, 'LineWidth', 2);
end

xlim([0.85 2.65]);
grid on; box on;
ylabel('DTWD','Interpreter','LaTex', 'FontSize',16);
set(gca,'xtick',[1:0.25:2.5],'xticklabel',names_)

end