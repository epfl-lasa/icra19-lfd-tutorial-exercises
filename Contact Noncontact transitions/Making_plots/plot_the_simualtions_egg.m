function plot_the_simualtions_egg(X,Center,Radiusy,Target,Option)
egg2=Option.egg2;
close all

%%
screensize = get( 0, 'Screensize' );
figure1 = figure();
set(figure1,'Position',screensize)
limits=Option.limits;
axes1 = axes('Parent',figure1);
% Create axes
hold(axes1,'on');
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
axis equal
xlabel('X [m]','Interpreter','latex');
ylabel('Y [m]','Interpreter','latex');
box(axes1,'on');
grid(axes1,'on');
set(axes1,'FontSize',20,'TickLabelInterpreter','latex');
ylim(axes1,[limits(3) limits(4)]);
xlim(axes1,[limits(1) limits(2)]);
th = linspace(0,2*pi) ;
if (egg2==1)
    y=Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    x= Radiusy*(sin(th))+Center(2);
else
    x=Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    y= Radiusy*(sin(th))+Center(2);
end
patch('YData',y,'XData',x,'FaceAlpha',0.6,...
    'LineStyle','none',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],'DisplayName','Contact surface') ;


plot(Center(1),Center(2),'DisplayName','Center of circle',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',20,...
    'LineStyle','none',...
    'Marker','o');


plot(Target(1),Target(2),'DisplayName','$x^*$',...
    'MarkerFaceColor',[0.494117647409439 0.184313729405403 0.556862771511078],...
    'MarkerEdgeColor','none',...
    'MarkerSize',20,...
    'LineStyle','none',...
    'Marker','hexagram');

if (egg2==1)
    t = linspace(0,2*pi);
    yin = Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    yout = (1+Option.rho)*((Radiusy)*(3+2*cos(th)-cos(th).*cos(th))/4-Radiusy/2)+Center(1);
    xin = Radiusy*(sin(th))+Center(2);
    xout= (Radiusy+Option.rho)*(sin(th))+Center(2);
else
    t = linspace(0,2*pi);
    xin = Radiusy*(3+2*cos(th)-cos(th).*cos(th))/4+Center(1)-Radiusy/2;
    xout = (1+Option.rho)*((Radiusy)*(3+2*cos(th)-cos(th).*cos(th))/4-Radiusy/2)+Center(1);
    yin = Radiusy*(sin(th))+Center(2);
    yout = (Radiusy+Option.rho)*(sin(th))+Center(2);
end
patch([xout,xin],[yout,yin],'g','linestyle','none','facealpha',0.3,'DisplayName','Transition region',...
    'FaceAlpha',0.3,...
    'LineStyle','none',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);


i=1;
plot(smooth(X{i}(1,:)),smooth(X{i}(2,:)),'LineWidth',2,'Color',[0 0.447058826684952 0.74117648601532],'DisplayName','Executed motion');
plot(smooth(X{i}(1,1)),smooth(X{i}(2,1)),'DisplayName','Initial Position',...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',5,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0.447058826684952 0.74117648601532]);
legend1 = legend(axes1,'show');
set(legend1,'Interpreter','latex');

for i=1:size(X,2)
    plot(smooth(X{i}(1,:)),smooth(X{i}(2,:)),'LineWidth',2,'Color',[0 0.447058826684952 0.74117648601532]);
    plot(smooth(X{i}(1,1)),smooth(X{i}(2,1)),'DisplayName','Initial Position',...
        'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
        'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
        'MarkerSize',5,...
        'Marker','o',...
        'LineStyle','none',...
        'Color',[0 0.447058826684952 0.74117648601532]);
end