function [A,Dem_i,check]=Construct_the_dynamcial_system(Poly,X_target,X_free,option)

clc
% close all
limits=option.limits;
X(1,:) = linspace(limits(1,1),limits(1,2),10^5);
X(2,:) = Poly(1)*X+Poly(2);
% screensize = get( 0, 'Screensize' );
% fig = figure();
% set(fig,'Position',screensize)
% plot(X(1,:),X(2,:),'DisplayName','The contact surface','LineWidth',4,...
%     'LineStyle','--',...
%     'Color',[0 0 0]);
hold on
% xlabel('X(1)','Interpreter','latex');
% ylabel('X(2)','Interpreter','latex');
plot(X_target(1,1),X_target(2,1),'DisplayName','The target','MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor','none',...
    'MarkerSize',30,...
    'Marker','pentagram',...
    'LineWidth',5,...
    'LineStyle','none');
legend1 = legend('show');
set(legend1,'Interpreter','latex','FontSize',20);
title('');
 legend('off') 
[Data,Dem_i]=generate_mouse_data(limits,1,1,'Draw some trajectories.');
legend('show')
for i=1:size(Dem_i,2)
    if sign(Poly(1)*X_free(1)+Poly(2)-X_free(2))*(Poly(1)*Dem_i(1,i)+Poly(2)-Dem_i(2,i))<0
        t = text(limits(1,1),(limits(1,4)+limits(1,3))/2,'ERROR: the demonstrations are started from the inside of the wall');
        s = t.FontSize;
        t.FontSize = 20;
        error('Program exit')
    end
end
d=2;
C=[];

Omega1=sdpvar(1,1);
A = [Omega1 0;0 Omega1];
C=C+[Omega1<= -0.0001];
Fun=(A*Data(1:d,:));
diff=Fun(1:d,:)-Data(d+1:2*d,:);

aux = sdpvar(d,length(diff));
Fun=sum((sum(aux.^2)));
C=C+[aux == diff];
options_solver=sdpsettings('solver','sedumi', ...
    'verbose', 0);
sol =  optimize(C,Fun,options_solver);
if sol.problem~=0
    disp('The optimization did not work. You need to change the solver or start it over!')
    check=0;
else
    disp('The dynamical system is successfully constructed!')
    check=1;
end
A= 10*value(A);
A = [zeros(d,d)  eye(d,d); A [-2*sqrt(-A(1,1)) 0;0 -2*sqrt(-A(2,2))]];