x1 = normrnd(5,2,10000,1); % First data set of 10000 values with mean 5
x2 = normrnd(10,2,10000,1); % Second data set of 10000 values with mean 10
x3 = normrnd(15,2,10000,1); % Third data set of 10000 values with mean 15

x = cat(2,x1,x2,x3); % Concatenate the data sets in a 10000 x 3 matrix

figure('Color',[1 1 1])
aboxplot(x,'labels',[5,10,15]); % Advanced box plot

xlabel('$\mu$','Interpreter','LaTex'); % Set the X-axis label

%% Grouped
% First group of 3 data set with standard deviation 2
x1 = normrnd(5,2,10000,1);
x2 = normrnd(10,2,10000,1);
x3 = normrnd(15,2,10000,1);
% Second group of 3 data set with standard deviation 4
y1 = normrnd(5,4,10000,1);
y2 = normrnd(10,4,10000,1);
y3 = normrnd(15,4,10000,1);
% Third group of 3 data set with standard deviation 6
z1 = normrnd(5,6,10000,1);
z2 = normrnd(10,6,10000,1);
z3 = normrnd(15,6,10000,1);

% Concatenate the data sets from each group in a 10000 x 3 matrix
x = cat(2,x1,x2,x3); 
y = cat(2,y1,y2,y3);
z = cat(2,z1,z2,z3);

size(x)
size(y)
size(z)
h = {x;y;z}; % Create a cell array with the data for each group

figure('Color',[1 1 1])
aboxplot(h,'labels',[5,10,15],'colorgrad','orange_down'); % Advanced box plot
grid on;
box on;
legend({'$\sigma=2$','$\sigma=4$','$\sigma=6$'},'Interpreter','LaTex')

xlabel('$\mu$','Interpreter','LaTex'); % Set the X-axis label

%%

  bar_input=rand(2,4)/2+0.5;
  errorbar_input=rand(2,4)/8;
  errorbar_groups(bar_input,errorbar_input, ...
      'bar_width',0.75,'errorbar_width',0.5, ...
      'optional_bar_arguments',{'LineWidth',1.5}, ...
      'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',1.5});
  
%%
A = [16 20 15 17 22 19 17]';
B = [22 15 16 16 16 18]';
C = [23 9 15 18 13 27 17 14 16 15 21 19 17]';

group = [    ones(size(A));
         2 * ones(size(B));
         3 * ones(size(C))];

figure
boxplot([A; B; C],group)
set(gca,'XTickLabel',{'A','B','C'})

  