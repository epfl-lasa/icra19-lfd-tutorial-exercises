function [Center,Radiusx,Radiusy,Target]=Construct_the_ellipsoid(options)

screensize = get( 0, 'Screensize' );
figure1 = figure();
set(figure1,'Position',screensize)
limits=options.limits;
axes1 = axes('Parent',figure1);
% Create axes
hold(axes1,'on');
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
axis equal
box(axes1,'on');
grid(axes1,'on');
set(axes1,'FontSize',20,'TickLabelInterpreter','latex');
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
title('Pick the center');
Center= ginput(1)';
disp('Specify the radius of the circle alog x axis.')
fig=scatter(Center(1),Center(2),150,'DisplayName','Center1','MarkerFaceAlpha',0.9,...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0]);
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
title('Specify the radius  x axis by adding one data point to that side.');
hold on
Radiusx= norm(ginput(1)'-Center);
disp('Specify the radius of the circle alog y axis.')
title('Specify the radius  y axis by adding one data point to that side.');
hold on
Radiusy= norm(ginput(1)'-Center);
% viscircles(X,R)
th = linspace(0,2*pi) ;
x = Radiusx*cos(th)+Center(1) ;
y = Radiusy*sin(th)+Center(2) ;
patch('YData',y,'XData',x,'FaceAlpha',0.6,...
    'LineStyle','none',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],'DisplayName','Contact surface') ;
delete(fig)
scatter(Center(1),Center(2),150,'DisplayName','Center','MarkerFaceAlpha',0.9,...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0]);

% axis equal
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
check=1;
while (check==1)
    title('Pick the target point as well, make sure that it is inside the surface');
    Target= ginput(1)';
    if (((Center(1)-Target(1))^2)/(Radiusx^2)+((Center(2)-Target(2))^2)/(Radiusy^2)-1<0)
        check=0;
    end
end
scatter(Target(1),Target(2),150,'DisplayName','Target','MarkerFaceAlpha',0.9,...
    'MarkerFaceColor',[0.494117647409439 0.184313729405403 0.556862771511078],...
    'MarkerEdgeColor','none',...
    'Marker','hexagram');

t = linspace(0,2*pi);
xin = Center(1) + Radiusx*cos(t);
xout =Center(1) + (Radiusx+options.rho)*cos(t);
yin = Center(2)+ Radiusy*sin(t);
yout = Center(2) + (Radiusy+options.rho)*sin(t);
patch([xout,xin],[yout,yin],'g','linestyle','none','facealpha',0.3,'DisplayName','Transition region',...
    'FaceAlpha',0.3,...
    'LineStyle','none',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);


disp('The contact surface is successfully constructed.')
legend(axes1,'show');
