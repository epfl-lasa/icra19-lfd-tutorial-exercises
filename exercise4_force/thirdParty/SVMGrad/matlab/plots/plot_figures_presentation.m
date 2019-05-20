%% Plot the training data
figure('Color',[1 1 1])
pos = find(labels == 1);
neg = find(labels == -1);
plot(X(pos,1), X(pos,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8); hold on;
plot(X(neg,1), X(neg,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8); hold on;

% Axis Labels
xlabel('$q_1^1$ (Joint 1-Robot 1)', 'Interpreter','Latex','FontSize',30); 
ylabel('$q_1^2$ (Joint 1 - Robot 2)', 'Interpreter','Latex','FontSize',30);
xt = get(gca, 'XTick');
set(gca,'FontSize', 30)
legend({'Non-Collided Joint Configurations', 'Collided Joint Configurations'}, 'Interpreter','Latex');

[hx,hy] = format_ticks(gca,{'$-\pi$','$-\frac{2\pi}{3}$','$-\frac{\pi}{3}$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$'},...
                           {'$-\pi$','$-\frac{\pi}{2}$','$\frac{\pi}{2}$','$\pi$'},...
                            [0,0.2,0.4,0.6,0.8,1], [0.4,0.6,0.8,1]);

%% Plot the training data + Boundary w/out gradient
limits = [0 1; 0.4 1];
plot_svmgrad_boundary(X, labels, svmgrad,  'draw', '', limits);
[hx,hy] = format_ticks(gca,{'$-\pi$','$-2/3\pi$','$-1/3\pi$','$1/3\pi$','$2/3\pi$','$\pi$'},...
                           {'$-\pi$','$-\frac{\pi}{2}$','$\frac{\pi}{2}$','$\pi$'},...
                            [0,0.2,0.4,0.6,0.8,1], [0.4,0.6,0.8,1]);
                        
                        
%% Plot the training data + Boundary w/out gradient
limits = [0 1; 0.4 1];
plot_svmgrad_boundary(X, labels, svmgrad,  'draw', 'grad', limits);
[hx,hy] = format_ticks(gca,{'$-\pi$','$-2/3\pi$','$-1/3\pi$','$1/3\pi$','$2/3\pi$','$\pi$'},...
                           {'$-\pi$','$-\frac{\pi}{2}$','$\frac{\pi}{2}$','$\pi$'},...
                            [0,0.2,0.4,0.6,0.8,1], [0.4,0.6,0.8,1]);