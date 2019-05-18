function [Center,Radiusy,Target]=Construct_the_egg(options)
egg2=options.egg2;
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
fig=scatter(Center(1),Center(2),150,'DisplayName','Center1','MarkerFaceAlpha',0.9,...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0]);
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
% viscircles(X,R)
th = linspace(0,2*pi) ;
if (egg2==1)
    disp('Specify the radius of the circle alog x axis.')
    title('Specify  the radius  x axis by adding one data point to that side.');
    hold on
    Radiusy= norm(ginput(1)'-Center);
    y=Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    x= Radiusy*(sin(th))+Center(2);
else
    disp('Specify the radius of the circle alog y axis.')
    title('Specify the radius  y axis by adding one data point to that side.');
    hold on
    Radiusy= norm(ginput(1)'-Center);
    
    x=Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    y= Radiusy*(sin(th))+Center(2);
end
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
    if (egg2==1)
        X_tmp=4*((Target(2)-Center(1))/(1+options.rho)+Radiusy/2)/Radiusy;
        Y_tmp=(Target(1)-Center(2))/(options.rho+Radiusy);
    else
        X_tmp=4*((Target(1)-Center(1))/(1+options.rho)+Radiusy/2)/Radiusy;
        Y_tmp=(Target(2)-Center(2))/(options.rho+Radiusy);
    end
    if ((X_tmp-3-Y_tmp^2)^2-4*(1-Y_tmp^2)<0)
        check=0;
    end
end
scatter(Target(1),Target(2),150,'DisplayName','Target','MarkerFaceAlpha',0.9,...
    'MarkerFaceColor',[0.494117647409439 0.184313729405403 0.556862771511078],...
    'MarkerEdgeColor','none',...
    'Marker','hexagram');


if (egg2==1)
    t = linspace(0,2*pi);
    yin = Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    yout = (1+options.rho)*((Radiusy)*(3+2*cos(th)-cos(th).*cos(th))/4-Radiusy/2)+Center(1);
    xin = Radiusy*(sin(th))+Center(2);
    xout= (Radiusy+options.rho)*(sin(th))+Center(2);
else
    t = linspace(0,2*pi);
    xin = Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    xout = (1+options.rho)*((Radiusy)*(3+2*cos(th)-cos(th).*cos(th))/4-Radiusy/2)+Center(1);
    yin = Radiusy*(sin(th))+Center(2);
    yout = (Radiusy+options.rho)*(sin(th))+Center(2);
end
patch([xout,xin],[yout,yin],'g','linestyle','none','facealpha',0.3,'DisplayName','Transition region',...
    'FaceAlpha',0.3,...
    'LineStyle','none',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);


disp('The contact surface is successfully constructed.')
legend(axes1,'show');
