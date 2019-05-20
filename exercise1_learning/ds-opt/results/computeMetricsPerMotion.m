function [] =  computeMetricsPerMotion(seds_rmse, seds_edot, seds_dtwd, em_lpv_rmse, em_lpv_edot, em_lpv_dtwd, pc_lpv_rmse, pc_lpv_edot, pc_lpv_dtwd,  names, data_type, face_color)


M = length(names);

for m=1:M
    fprintf('**************** Performance of Different Algorithms on (%s)-Motion ****************\n', names{m});
    seds_rmse_mean = mean(seds_rmse(:,m));
    seds_rmse_std  = std(seds_rmse(:,m));
    fprintf('SEDS RMSE on %s Set %2.4f +/- %2.4f\n', data_type, seds_rmse_mean, seds_rmse_std);
    
    seds_edot_mean = mean(seds_edot(:,m));
    seds_edot_std  = std(seds_edot(:,m));
    fprintf('SEDS e-dot on %s Set %2.4f +/- %2.4f\n', data_type, seds_edot_mean, seds_edot_std);
    
    seds_dtwd_mean = mean(seds_dtwd(:,m));
    seds_dtwd_std  = std(seds_dtwd(:,m));
    fprintf('SEDS DTWD on %s Set %2.4f +/- %2.4f\n', data_type, seds_dtwd_mean, seds_dtwd_std);
    
    opt_em = size(em_lpv_rmse,3);
    for opt=1:opt_em
        em_lpv_rmse_mean = mean(em_lpv_rmse(:,m,opt));
        em_lpv_rmse_std  = std(em_lpv_rmse(:,m,opt));
        fprintf('EM-GMM LPV (O%d) RMSE on %s Set %2.4f +/- %2.4f\n', opt, data_type, em_lpv_rmse_mean, em_lpv_rmse_std);
        
        em_lpv_edot_mean = mean(em_lpv_edot(:,m,opt));
        em_lpv_edot_std  = std(em_lpv_edot(:,m,opt));
        fprintf('EM-GMM LPV (O%d) e-dot on %s Set %2.4f +/- %2.4f\n', opt, data_type, em_lpv_edot_mean, em_lpv_edot_std);
        
        em_lpv_dtwd_mean = mean(em_lpv_dtwd(:,m,opt));
        em_lpv_dtwd_std  = std(em_lpv_dtwd(:,m,opt));
        fprintf('EM-GMM LPV (O%d) DTWD on %s Set %2.4f +/- %2.4f\n', opt, data_type, em_lpv_dtwd_mean, em_lpv_dtwd_std);        
    end
    
    opt_pc = size(pc_lpv_rmse,3);
    for opt=1:opt_pc
        pc_lpv_rmse_mean = mean(pc_lpv_rmse(:,m,opt));
        pc_lpv_rmse_std  = std(pc_lpv_rmse(:,m,opt));
        fprintf('PC-GMM LPV (O%d) RMSE on %s Set %2.4f +/- %2.4f\n', opt, data_type, pc_lpv_rmse_mean, pc_lpv_rmse_std);
        
        pc_lpv_edot_mean = mean(pc_lpv_edot(:,m,opt));
        pc_lpv_edot_std  = std(pc_lpv_edot(:,m,opt));
        fprintf('PC-GMM LPV (O%d) e-dot on %s Set %2.4f +/- %2.4f\n', opt, data_type, pc_lpv_edot_mean, pc_lpv_edot_std);
        
        pc_lpv_dtwd_mean = mean(pc_lpv_dtwd(:,m,opt));
        pc_lpv_dtwd_std  = std(pc_lpv_dtwd(:,m,opt));
        fprintf('PC-GMM LPV (O%d) DTWD on %s Set %2.4f +/- %2.4f\n', opt, data_type, pc_lpv_dtwd_mean, pc_lpv_dtwd_std);
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
%      edge_C = [rgb('MidNightBlue');
%          rgb('Indigo');
%          rgb('DarkSlateBlue');
%          rgb('BlueViolet');
%          rgb('RoyalBlue');
%          rgb('DodgerBlue');
%          rgb('SteelBlue');
%          ];
%      
%      face_C = [rgb('MediumBlue');
%          rgb('DarkMagenta');
%          rgb('DarkOrchid');
%          rgb('DarkViolet');
%          rgb('DeepSkyBlue');
%          rgb('SkyBlue');
%          rgb('LightSkyBlue');
%          ];
%      
    Y_rmse = [seds_rmse_mean mean(em_lpv_rmse(:,m,1)) mean(em_lpv_rmse(:,m,2)) mean(em_lpv_rmse(:,m,3)) ...
        mean(pc_lpv_rmse(:,m,1)) mean(pc_lpv_rmse(:,m,2)) mean(pc_lpv_rmse(:,m,3))];
    E_rmse = 2*[seds_rmse_std std(em_lpv_rmse(:,m,1)) std(em_lpv_rmse(:,m,2)) std(em_lpv_rmse(:,m,3)) ...
        std(pc_lpv_rmse(:,m,1)) std(pc_lpv_rmse(:,m,2)) std(pc_lpv_rmse(:,m,3))];
    X = [1 1.25 1.5 1.75 2 2.25 2.5];   
    
    if face_color
        h = superbar(X, Y_rmse, 'BarFaceColor', face_C,'BarEdgeColor', edge_C, 'E', E_rmse, 'ErrorbarStyle', 'T');
        set(h, 'LineWidth', 1);
    else 
        h = superbar(X, Y_rmse, 'BarFaceColor', 'none','BarEdgeColor', edge_C, 'E', E_rmse, 'ErrorbarStyle', 'T');
        set(h, 'LineWidth', 2);
    end
    
    xlim([0.85 2.65]);    
    names_ = {'S','E(O1)','E(O2)','E(O3)', 'PC(O1)', 'PC(O2)', 'PC(O3)'};
    set(gca,'xtick',[1:0.25:2.5],'xticklabel',names_)
    grid on; box on;
    ylabel('RMSE','Interpreter','LaTex', 'FontSize',20)
    title(sprintf('Performance on $%s$ %s-set',names{m}, data_type), 'Interpreter','LaTex', 'FontSize',18)
    

    subplot(3,1,2);
    Y_edot = [seds_edot_mean mean(em_lpv_edot(:,m,1)) mean(em_lpv_edot(:,m,2)) mean(em_lpv_edot(:,m,3)) ...
        mean(pc_lpv_edot(:,m,1)) mean(pc_lpv_edot(:,m,2)) mean(pc_lpv_edot(:,m,3))];
    E_edot = 2*[seds_edot_std std(em_lpv_edot(:,m,1)) std(em_lpv_edot(:,m,2)) std(em_lpv_edot(:,m,3)) ...
        std(pc_lpv_edot(:,m,1)) std(pc_lpv_edot(:,m,2)) std(pc_lpv_edot(:,m,3))];
    if face_color
        h = superbar(X, Y_edot, 'BarFaceColor', face_C,'BarEdgeColor', edge_C,'E', E_edot, 'ErrorbarStyle', 'T');
        set(h, 'LineWidth', 1);
    else
        h = superbar(X, Y_edot, 'BarFaceColor', 'none','BarEdgeColor', edge_C, 'E', E_edot, 'ErrorbarStyle', 'T');
        set(h, 'LineWidth', 2);
    end
    
    xlim([0.85 2.65]);    
    grid on; box on;
    ylabel('$\dot{e}$','Interpreter','LaTex', 'FontSize',20)
    set(gca,'xtick',[1:0.25:2.5],'xticklabel',names_)
    
    Y_dtwd = [seds_dtwd_mean mean(em_lpv_dtwd(:,m,1)) mean(em_lpv_dtwd(:,m,2)) mean(em_lpv_dtwd(:,m,3)) ...
        mean(pc_lpv_dtwd(:,m,1)) mean(pc_lpv_dtwd(:,m,2)) mean(pc_lpv_dtwd(:,m,3))];
    E_dtwd = 2*[seds_dtwd_std std(em_lpv_dtwd(:,m,1)) std(em_lpv_dtwd(:,m,2)) std(em_lpv_dtwd(:,m,3)) ...
        std(pc_lpv_dtwd(:,m,1)) std(pc_lpv_dtwd(:,m,2)) std(pc_lpv_dtwd(:,m,3))];    
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
    ylabel('DTWD','Interpreter','LaTex', 'FontSize',20);           
    set(gca,'xtick',[1:0.25:2.5],'xticklabel',names_)
end

%     violinplot([seds_rmse(:,m) em_lpv_rmse(:,m,1) em_lpv_rmse(:,m,2) em_lpv_rmse(:,m,3) pc_lpv_rmse(:,m,1) pc_lpv_rmse(:,m,2) pc_lpv_rmse(:,m,3)], ...
%          {'S','E(O1)','E(O2)','E(O3)', 'PC(O1)', 'PC(O2)', 'PC(O3)'});          
%     violinplot([seds_edot(:,m) em_lpv_edot(:,m,1) em_lpv_edot(:,m,2) em_lpv_edot(:,m,3) pc_lpv_edot(:,m,1) pc_lpv_edot(:,m,2) pc_lpv_edot(:,m,3)], ...
%         {'S','E(O1)','E(O2)','E(O3)', 'PC(O1)', 'PC(O2)', 'PC(O3)'});
%     violinplot([seds_dtwd(:,m) em_lpv_dtwd(:,m,1) em_lpv_dtwd(:,m,2) em_lpv_dtwd(:,m,3) pc_lpv_dtwd(:,m,1) pc_lpv_dtwd(:,m,2) pc_lpv_dtwd(:,m,3)], ...
%         {'S','E(O1)','E(O2)','E(O3)', 'PC(O1)', 'PC(O2)', 'PC(O3)'});

end