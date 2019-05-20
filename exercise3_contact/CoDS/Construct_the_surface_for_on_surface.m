function [N_x,p,X_Target_projections,X_free,check]=Construct_the_surface_for_on_surface(options)

% screensize = get( 0, 'Screensize' );
%  fig = figure();
%  axes(options.Figure);
%  set(fig,'Position',screensize)
limits=options.limits;
% fig=options.Figure;
%  limits=[limits(1) limits(2) limits(3) limits(4)];
[Datat_Wall,~]=generate_mouse_data(limits,0,0,'Draw the wall');
p = polyfit(Datat_Wall(1,:),Datat_Wall(2,:),1);
X(1,:) = linspace(limits(1,1),limits(1,2),10^5);
X(2,:) = p(1)*X+p(2);
% close all
% pause(0.5)
disp('Specify the free motion side.')
delete(findobj(gca, 'type', 'line'));
createfigure_with_wall(X(1,:),X(2,:),'Specify the free motion side by adding one data point to that side.');
X_free= ginput(1)';
if abs(p(1)*X_free(1)+p(2)-X_free(2))<0.1
    t = text(limits(1,1),(limits(1,4)+limits(1,3))/2,{'ERROR: the chosen point is very close to the contact surface.';'Please pick one point which is far from the contact surface.'});
    s = t.FontSize;
    t.FontSize = 20;
    X_free= ginput(1)';
    if abs(p(1)*X_free(1)+p(2)-X_free(2))<0.1
        close all
        clc
        t = text(limits(1,1),(limits(1,4)+limits(1,3))/2,{'ERROR: the chosen point is very close to the contact surface.';'Please start again and pick one point which is far from the contact surface.'});
        s = t.FontSize;
        t.FontSize = 20;
%         pause(1)
        error('ERROR: the chosen point is very close to the contact surface. Please start again and pick one point which is far from the contact surface.')
    end
end
% X_free=[0 ;0];
[~,Distance_i]=min(sum((repmat(X_free(:,1),1,size(X,2))-X).*...
    (repmat(X_free(:,1),1,size(X,2))-X),1));

handle=sign(X_free-X(:,Distance_i));

N_x=[handle(1,1)*abs(p(1));handle(2,1)*1]/norm([p(1);1]);

% close all
% pause(0.5)
disp('Specify the target position.')
createfigure_with_wall(X(1,:),X(2,:),'Specify the target position by adding one data point on the contact surface.');
X_Contact= ginput(1)';
[~,Distance_i]=min(sum((repmat(X_Contact(:,1),1,size(X,2))-X).*...
    (repmat(X_Contact(:,1),1,size(X,2))-X),1));

X_Target_projections=X(:,Distance_i);

disp('The contact surface is successfully constructed.')
check=1;