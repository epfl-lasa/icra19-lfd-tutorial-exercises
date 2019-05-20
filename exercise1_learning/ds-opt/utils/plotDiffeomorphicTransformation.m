function [h_fig] = plotDiffeomorphicTransformation(nDemos, source, target_trajectory, alltarget_trajectory, allSource, diff_fun, inverse_diff_fun)
% Create figure
h_fig = figure('Color',[1 1 1]); hold all;

% Demonstrated Trajectory transformed to Virtual through \phi^-1(\xi)
invPhiTarg = inverse_diff_fun(target_trajectory);

% Virtual Trajectory transformed to Demonstration through \\phi(\chi)
PhiSrc = diff_fun(source);

dim = size(source,1);
if dim == 2
    %Max for plotting
    XLimPlot = [min([1.2*alltarget_trajectory, 0.8*alltarget_trajectory, -.25*ones(dim,1)],[], 2), max([1.2*alltarget_trajectory, 0.8*alltarget_trajectory, .25*ones(dim,1)],[],2)];
    XLimPlot2 = [min([1.2*allSource, 0.8*allSource, -.25*ones(dim,1)],[], 2), max([1.2*allSource, 0.8*allSource, .25*ones(dim,1)],[],2)];
    
    % Demonstrated Mean Trajectory \xi
    scatter(target_trajectory(1,:),target_trajectory(2,:), 20, [1 0 0],'filled');
    if nDemos > 1
        % All Demonstrated Trajectories \xi
        scatter(alltarget_trajectory(1,:),alltarget_trajectory(2,:), 10, [0 0 0],'filled');
    end
    % Virtual Trajectory \chi
    scatter(source(1,:),source(2,:), 20, [0 0 1], 'filled'); hold on;
    scatter(invPhiTarg(1,:), invPhiTarg(2,:), 20, [0 0.5 1], '+'); hold on;
    scatter(PhiSrc(1,:), PhiSrc(2,:), 20, [1 0.5 0], '+'); hold on;
    
    % Plotting the deformed grid through the diffeomorphism
    plotGrid(1, XLimPlot2, 10, 1000, diff_fun);
    xlim(XLimPlot(1,:));
    ylim(XLimPlot(2,:));
    xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
    ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
    
elseif dim == 3
    
    % Demonstrated Mean Trajectory \xi
    scatter3(target_trajectory(1,:),target_trajectory(2,:),target_trajectory(3,:), 20, [1 0 0],'filled');
    if nDemos > 1
        % All Demonstrated Trajectories \xi
        scatter3(alltarget_trajectory(1,:),alltarget_trajectory(2,:),alltarget_trajectory(3,:), 10, [0 0 0],'filled');
    end
    % Virtual Trajectory \chi
    scatter3(source(1,:),source(2,:), source(3,:), 20, [0 0 1], 'filled'); hold on;
    scatter3(invPhiTarg(1,:), invPhiTarg(2,:), invPhiTarg(3,:), 20, [0 0.5 1], '+'); hold on;
    scatter3(PhiSrc(1,:), PhiSrc(2,:),PhiSrc(3,:), 20, [1 0.5 0], '+'); hold on;
    xlabel('$\xi_1$','Interpreter','LaTex','FontSize',20);
    ylabel('$\xi_2$','Interpreter','LaTex','FontSize',20);
    zlabel('$\xi_3$','Interpreter','LaTex','FontSize',20);
 
end
grid on;
if nDemos > 1
    legend({'Mean Demo trajectory $\Xi=\{\xi_1,\dots,\xi_T\}$','Demonstrated trajectories $\Xi=\{\xi_1,\dots,\xi_T\}$', ... ,
        'Virtual Linear Trajectory $\chi=\{\chi_1,\dots,\chi_T\}$', ...
        'Transformed Demonstrated Trajectory $\Phi^{-1}(\Xi)=\{\phi^{-1}(\xi_1),\dots,\phi^{-1}(\xi_T)\}$',...
        'Transformed Virtual Trajectory $\Phi(\chi)=\{\phi(\chi_1),\dots,\phi(\chi_T)\}$'},'Interpreter','LaTex','FontSize',10)
else
    legend({'Mean Demo trajectory $\Xi=\{\xi_1,\dots,\xi_T\}$', ... ,
        'Virtual Linear Trajectory $\chi=\{\chi_1,\dots,\chi_T\}$', ...
        'Transformed Demonstrated Trajectory $\Phi^{-1}(\Xi)=\{\phi^{-1}(\xi_1),\dots,\phi^{-1}(\xi_T)\}$',...
        'Transformed Virtual Trajectory $\Phi(\chi)=\{\phi(\chi_1),\dots,\phi(\chi_T)\}$'},'Interpreter','LaTex','FontSize',10)
end
title('Results of Diffeomorphic Matching Algorithm','Interpreter','LaTex', 'FontSize',15)
       
if dim == 3
    view(-120, 30);
    axis equal
end


end