function [h, vals]  = plot_svmgrad_boundary_lines( data, labels, model, grad_points)
% PLOT_SVM_BOUNDARY plots the training data
%   and decision boundary, given a model produced by LIBSVM
%   input ----------------------------------------------------------------
%
%       o data      : (1 x 1), number of data points to generate.
%
%       o labels    : (1 x 1), dimension of the data. [2 or 3]
%
%       o model     : (1 x 1), number of classes.
%
%       o options   : struct
%
%       o varargin  : string, if 'draw' draw contours otherwise, don't
%

% Check Labels
labels(find(labels==0)) = -1;

% Create Figure
figure('Color',[1 1 1]); 
hold on

% Make classification predictions over a grid of values
xplot = linspace(min(data(:,1)), max(data(:,1)), 100)';
yplot = linspace(min(data(:,2))-0.1, max(data(:,2)), 100)';
[X, Y] = meshgrid(xplot, yplot);
vals = zeros(size(X));

% Using SVMGrad Library
for i = 1:size(X, 2)   
   X_row = [X(:,i),Y(:,i)]'; 
   values    = zeros(length(X_row),1);
   for ii=1:length(X_row)
        x = X_row(:,ii)';
        values(ii,:)    = calculateGamma( model, x' );
   end
   vals(:,i) = values;
end

% Plot the SVM Contours
level = 20; n = ceil(level/2);
cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
cmap = [cmap1; cmap2(2:end, :)];
colormap(vivid(cmap, [.5, .5]));
contourf(X,Y, vals, 50, 'LineStyle', 'none');
colorbar


% Plot the SVM Decision Boundaries
contour(X,Y, vals, [0 0], 'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
contour(X,Y, vals, [1 1], 'LineWidth', 3, 'Color', 'k');
contour(X,Y, vals, [2 2], 'LineWidth', 2, 'Color', 'k');
contour(X,Y, vals, [3 3], 'LineWidth', 2, 'Color', 'k');
contour(X,Y, vals, [4 4], 'LineWidth', 2, 'Color', 'k');
contour(X,Y, vals, [5 5], 'LineWidth', 2, 'Color', 'k');
contour(X,Y, vals, [6 6], 'LineWidth', 2, 'Color', 'k');
contour(X,Y, vals, [7 7], 'LineWidth', 2, 'Color', 'k');
hold on;

% Plot the training data on top of the boundary
pos = find(labels == 1);
neg = find(labels == -1);
plot(data(pos,1), data(pos,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 5); hold on;
plot(data(neg,1), data(neg,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 5); hold on;

% Extract point for gradient computation
U = zeros(size(grad_points,1),1);
V = zeros(size(grad_points,1),1);
for i = 1:size(grad_points, 1)    
   gradient = calculateGammaDerivative( model, grad_points(i,:)' );
   gradient = gradient/norm(gradient);
   U(i,1)   = gradient(1);
   V(i,1)   = gradient(2);
end
quiver(grad_points(:,1), grad_points(:,2), U, V, 0.25,  'Color', 'k', 'LineWidth',4);
text(0.6231,0.7461,'6.967','BackgroundColor','w')

% Plot Points for Gradient Vectors
scatter(grad_points(:,1),grad_points(:,2),100,'d','MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 0], 'LineWidth', 2);

% Axis Labels
% xlabel('q_1^1 (Joint 1-Robot 1)', 'FontSize',30); ylabel('q_1^2 (Joint 1 - Robot 2)', 'FontSize',30);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)

% Legends
% legend_names = {'\Gamma(q^{ij}) Self-Collision Boundary Function','\Gamma(q^{ij}) = 0 (Classifier Boundary)', '\Gamma(q^{ij})= +1 (Self-Collision Boundary)', '\nabla\Gamma(q^{ij}) (Gradient of Boundary Function)'};
% legend(legend_names,'Location','NorthWest', 'FontSize',20);

axis equal
grid on
box on

end

