function plot_the_simualtions(DDX,DX,X,F,Time,Poly,X_initial,X_target,X_free,X_C,X_L,N_x,Option)

% close all

Limits=Option.limits;

if (Option.animation==1)
    plot_animation(X,Poly,X_initial,X_target,X_free,X_C,X_L,N_x,Option)
end
%%
% screensize = get( 0, 'Screensize' );
% fig = figure();
% subplot1 = subplot(1,1,1);
% set(fig,'Position',screensize)
cla
X1=[ Option.axes1.XLim(1):0.05: Option.axes1.XLim(2)];
X2=[ Option.axes1.YLim(1):0.05: Option.axes1.YLim(2)]';

X1=repmat(X1,size(X2,1),1);
X2=repmat(X2,1,size(X1,2));
Walla=zeros(size(X2));
Wall_Base=Option.N_x'*Option.X_C;

Handle_sign=sign(Option.N_x'*Option.X_free-Wall_Base);
for ii=1:size(X2,1),
    for jj=1:size(X2,2)
        XX=[X1(ii,jj);X2(ii,jj)];
        Walla(ii,jj)=Handle_sign*(Option.N_x'*XX-Wall_Base)+...
            (Option.rho-(Option.X_L-Option.X_C)'*(Option.X_L-XX))*exp(-Option.kamma_slider*(Option.X_L-XX)'*(Option.X_L-XX));
        if  Option.rho<(Walla(ii,jj))
            Walla(ii,jj)=Option.rho;
        end
        if Walla(ii,jj)<-2
            Walla(ii,jj)=-2;
        end
    end
end
clim=[-2 Option.rho];
hold on
contourf(X1,X2,Walla);
colormap(hot);

% plot_Wall_counters(subplot1,N_x,X_free,X_C,X_L,Option.rho,Limits)
h2=plot(X_target(1,1),X_target(2,1),'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','none',...
    'MarkerSize',30,...
    'Marker','pentagram',...
    'LineWidth',5,...
    'LineStyle','none');
hold on
h1=plot(X_initial(1,:),X_initial(2,:),...
    'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
    'MarkerEdgeColor','none',...
    'MarkerSize',30,...
    'Marker','hexagram',...
    'LineWidth',5,...
    'LineStyle','none');
X_wall(1,:) = linspace(Limits(1,1),Limits(1,2),10^5);
X_wall(2,:) = Poly(1)*X_wall+Poly(2);
h3=plot(X_wall(1,:),X_wall(2,:),'LineWidth',4,...
    'LineStyle','--',...
    'Color',[0 0 0]);
hold on
h4=plot(X_C(1,1),X_C(2,1),...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',30,...
    'Marker','^',...
    'LineStyle','none');
for i=1:size(X,2)
    h8{i} = plot(X{i}(1,:),X{i}(2,:),'LineWidth',4,'Color',[0 0 i/size(X,2)]);
end
if Option.Onsurface==0
    h5=plot(X_L(1,1),X_L(2,1),...
        'MarkerFaceColor',[1 0 0],...
        'MarkerSize',30,...
        'Marker','v',...
        'LineStyle','none');
end
grid on 
box on


% legend2=legend([h1(1) h2 h3 h4 h5 h8{i}],'The initial positions','The target','The contact surface', 'Desired contact location', 'Leaving point','Executed motion');
% set(legend2,'Orientation','vertical','Location','best');
%%
% clearvars h1 h2 h3 h4 h5 h6 h7 h8
% subplot1 = subplot(5,2,[8 10]);
% hold(subplot1,'on');
% box(subplot1,'on');
% grid(subplot1,'on');
% title('Velocity of the motion normal to the surface.')
% xlabel('$Time [s]$','Interpreter','latex');
% ylabel('$Proj_{\dot{X}}$','Interpreter','latex');
% set(subplot1,'FontSize',20);
% for i=1:size(X,2)
%     h1{i} = plot(Time{i}',N_x'*DX{i},'LineWidth',4,'Color',[0 0 i/size(X,2)]);
% end
%%
% clearvars h1 h2
% subplot1 = subplot(5,2,[7 9]);
% hold(subplot1,'on');
% box(subplot1,'on');
% grid(subplot1,'on');
% set(subplot1,'FontSize',20);
% title('Force exerted by the surface to the motion.')
% for i=1:size(X,2)
%     h1{i} = plot(Time{i}',N_x'*F{i},'LineWidth',4,'Color',[0 0 i/size(X,2)]);
% end
% xlabel('$Time [s]$','Interpreter','latex');
% ylabel('$F[N]$','Interpreter','latex');





function plot_animation(X,Poly,X_initial,X_target,X_free,X_C,X_L,N_x,Option)
Limits=Option.limits;
% close all;
% screensize = get( 0, 'Screensize' );
% fig = figure();
% axes1 = axes('Parent',fig);
% set(fig,'Position',screensize)
% plot_Wall_counters(axes1,N_x,X_free,X_C,X_L,Option.rho,Limits)
X1=[ Option.axes1.XLim(1):0.05: Option.axes1.XLim(2)];
X2=[ Option.axes1.YLim(1):0.05: Option.axes1.YLim(2)]';
cla
legend('off');
X1=repmat(X1,size(X2,1),1);
X2=repmat(X2,1,size(X1,2));
Walla=zeros(size(X2));
Wall_Base=Option.N_x'*Option.X_C;

Handle_sign=sign(Option.N_x'*Option.X_free-Wall_Base);
for ii=1:size(X2,1),
    for jj=1:size(X2,2)
        XX=[X1(ii,jj);X2(ii,jj)];
        Walla(ii,jj)=Handle_sign*(Option.N_x'*XX-Wall_Base)+...
            (Option.rho-(Option.X_L-Option.X_C)'*(Option.X_L-XX))*exp(-Option.kamma_slider*(Option.X_L-XX)'*(Option.X_L-XX));
        if  Option.rho<(Walla(ii,jj))
            Walla(ii,jj)=Option.rho;
        end
        if Walla(ii,jj)<-2
            Walla(ii,jj)=-2;
        end
    end
end
clim=[-2 Option.rho];
contourf(X1,X2,Walla);
colormap(hot);
hold on

h2=plot(X_target(1,1),X_target(2,1),'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','none',...
    'MarkerSize',30,...
    'Marker','pentagram',...
    'LineWidth',5,...
    'LineStyle','none');

h1=plot(X_initial(1,:),X_initial(2,:),...
    'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
    'MarkerEdgeColor','none',...
    'MarkerSize',30,...
    'Marker','hexagram',...
    'LineWidth',5,...
    'LineStyle','none');
hold on
X_wall(1,:) = linspace(Limits(1,1),Limits(1,2),10^5);
X_wall(2,:) = Poly(1)*X_wall+Poly(2);
h3=plot(X_wall(1,:),X_wall(2,:),'LineWidth',4,...
    'LineStyle','--',...
    'Color',[0 0 0]);
hold on

h4=plot(X_C(1,1),X_C(2,1),...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',30,...
    'Marker','^',...
    'LineStyle','none');
if Option.Onsurface==0
    h5=plot(X_L(1,1),X_L(2,1),...
        'MarkerFaceColor',[1 0 0],...
        'MarkerSize',30,...
        'Marker','v',...
        'LineStyle','none');
    legend2=legend([h1(1) h2 h3 h4 h5],'The initial positions','The target','The contact surface', 'Desired contact location', 'Leaving point');
else
    legend2=legend([h1(1) h2 h3 h4],'The initial positions','The target','The contact surface', 'Desired contact location');
end
set(legend2,'Orientation','horizontal','Location','northwest');
SIZE=size(X{1},2);
for i=1:size(X,2)
    H_animated{i}=animatedline('Color','b','LineWidth',3);
    if size(X{i},2)<SIZE
        SIZE=size(X{i},2);
    end
end
counter=1;
while counter<SIZE
    for i=1:size(X,2)
        addpoints(H_animated{i},X{i}(1,counter),X{i}(2,counter))
        drawnow
    end
    counter=int64(counter+size(X{i},2)/abs(10*Option.Tfinal));
end
