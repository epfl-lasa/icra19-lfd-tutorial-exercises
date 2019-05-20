function [h, vals]  = plot_svmgrad_boundary( data, labels, model, varargin)
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

if (size(varargin, 2) <= 2)
    % Default limits
    limits = [min(data(:,1)), max(data(:,1));min(data(:,2))-0.1, max(data(:,2))];
else
    % User-defined limits
    limits = varargin{3};
end

% Make classification predictions over a grid of values
xplot = linspace(limits(1,1), limits(1,2), 100)';
yplot = linspace(limits(2,1), limits(2,2), 100)';
[X, Y] = meshgrid(xplot, yplot);
vals = zeros(size(X));

U = zeros(size(X));
V = zeros(size(Y));

% Using SVMGrad Library
for i = 1:size(X, 2)   
   X_row = [X(:,i),Y(:,i)]'; 
   values    = zeros(length(X_row),1);
   gradients = zeros(length(X_row),2);   
   for ii=1:length(X_row)
        x = X_row(:,ii)';
        values(ii,:)    = calculateGamma( model, x' );
        gradients(ii,:) = calculateGammaDerivative( model, x' )';
   end
   vals(:,i) = values;
   U(:,i)    = gradients(:,1);
   V(:,i)    = gradients(:,2);
end

% Plot the SVM Contours
level = 20; n = ceil(level/2);
cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
cmap = [cmap1; cmap2(2:end, :)];
colormap(vivid(cmap, [.5, .5]));
if (size(varargin, 2) > 0) && strcmp(varargin{1}, 'draw')
    contourf(X,Y, vals, 50, 'LineStyle', 'none');   
%     colorbar
elseif (size(varargin, 2) > 0) && strcmp(varargin{1}, 'surf')
    surf(X,Y, vals); shading interp;
    colorbar
else       
end

% Plot the SVM Decision Boundaries
% contour(X,Y, vals, [0 0], 'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
contour(X,Y, vals, [1 1], 'LineWidth', 3, 'Color', 'k');
hold on;

% Extract some learnt parameters and options
SVs              = model.SVs;

% Plot Support Vectors
% scatter(SVs(1,:),SVs(2,:),30,'o','MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [1 1 0], 'LineWidth', 2);

% Extract point for gradient computation
if (size(varargin, 2) == 2) && strcmp(varargin{2}, 'grad')
    [id] = find(vals>1);
    rand_id = randperm(length(id));
    grad_id = id(rand_id(1:round(length(id)/3)));
    quiver(X(grad_id), Y(grad_id), U(grad_id), V(grad_id), 'Color', 'k', 'AutoScale','on', 'LineWidth',2);
end

% Plot the training data on top of the boundary
pos = find(labels == 1);
neg = find(labels == -1);
plot(data(pos,1), data(pos,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 6); hold on;
plot(data(neg,1), data(neg,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;

% Axis Labels
xlabel('$q_1^1$ (Joint 1-Robot 1)', 'Interpreter','Latex','FontSize',30); 
ylabel('$q_1^2$ (Joint 1 - Robot 2)', 'Interpreter','Latex','FontSize',30);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)

% Legends
% legend_names = {'\Gamma(q^{ij}) Self-Collision Boundary Function','\Gamma(q^{ij}) = 0 (Classifier Boundary)', '\Gamma(q^{ij})= +1 (Self-Collision Boundary)','y = +1 (Free Robot Configurations)','y = -1 (Collided Robot Configurations)', '\nabla\Gamma(q^{ij}) (Gradient of Boundary Function)'};
% legend_names = {'\Gamma(q^{ij}) = 0 (Classifier Boundary)', '\Gamma(q^{ij})= +1 (Self-Collision Boundary)','y = +1 (Free Robot Configurations)','y = -1 (Collided Robot Configurations)', '\nabla\Gamma(q^{ij}) (Gradient of Boundary Function)'};
% legend_names = {'\Gamma(q^{ij}) Self-Collision Boundary Function','\Gamma(q^{ij}) = 0 (Classifier Boundary)', '\Gamma(q^{ij})= +1 (Self-Collision Boundary)', '\nabla\Gamma(q^{ij}) (Gradient of Boundary Function)'};
% legend(legend_names,'Location','NorthWest', 'FontSize',20);

axis equal
grid on
box on

end

