function [N_x,p,X_target,check]=Construct_the_surface(options)

% screensize = get( 0, 'Screensize' );
% fig = figure();
% set(fig,'Position',screensize)
limits=options.limits;
% limits=[limits(1) limits(2) 0 limits(4)];
 legend('off');
[Datat_Wall,~]=generate_mouse_data(limits,0,0,'Draw the wall');
 legend('show');
p = polyfit(Datat_Wall(1,:),Datat_Wall(2,:),1);
X(1,:) = linspace(limits(1,1),limits(1,2),10^5);
X(2,:) = p(1)*X+p(2);
% close all
delete(findobj(gca, 'type', 'line'));
% pause(0.5)
disp('Specify the target position.')
createfigure_with_wall(X(1,:),X(2,:),'Specify the target point by adding one data point to that side.')
X_free= ginput(1)';
X_target=X_free;
% X_free=[0 ;0];
[~,Distance_i]=min(sum((repmat(X_free(:,1),1,size(X,2))-X).*...
    (repmat(X_free(:,1),1,size(X,2))-X),1));

handle=sign(X_free-X(:,Distance_i));

N_x=[handle(1,1)*abs(p(1));handle(2,1)*1]/norm([p(1);1]);

if ((1-N_x'*(X_free-X(:,Distance_i))/norm(X_free-X(:,Distance_i)))<0.0001)
    disp('The contact surface is successfully constructed.')
    check=1;
else
    disp('Something is wrong, please start it over or contact Sina Mirrazavi')
    check=0;
end