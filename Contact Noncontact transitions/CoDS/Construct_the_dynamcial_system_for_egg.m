function [A,Dem_i,check]=Construct_the_dynamcial_system_for_egg(Center,Radiusy,Option)
egg2=Option.egg2;
clc
limits=Option.limits;
hold on
% xlabel('X(1)','Interpreter','latex');
% ylabel('X(2)','Interpreter','latex');
legend1 = legend('off');
set(legend1,'Interpreter','latex','FontSize',20);
[Data,Dem_i]=generate_mouse_data(limits,1,1,'Draw some motions strating from the free-space region, make sure that it ends up at the target point and it goes through the contact surface !');

for i=1:size(Dem_i,2)
    if (egg2==1)
        X_tmp=4*((Dem_i(2,i)-Center(1))/(1+Option.rho)+Radiusy/2)/Radiusy;
        Y_tmp=(Dem_i(1,i)-Center(2))/(Option.rho+Radiusy);
    else
        X_tmp=4*((Dem_i(2,i)-Center(1))/(1+Option.rho)+Radiusy/2)/Radiusy;
        Y_tmp=(Dem_i(1,i)-Center(2))/(Option.rho+Radiusy);
    end
    if ((X_tmp-3-Y_tmp^2)^2-4*(1-Y_tmp^2)<0)
        t = text(limits(1,1),(limits(1,4)+limits(1,3))/2,'ERROR: the demonstrations are started from the inside of the wall');
        s = t.FontSize;
        t.FontSize = 20;
        error('Program exit')
    end
end
d=2;
C=[];

Omega1=sdpvar(1,1);
Omega2=sdpvar(1,1);
A = [Omega1 0;0 Omega2];
C=C+[Omega1<= -0.0001]+[Omega2<= -0.0001];
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
A=10* value(A);
A = [zeros(d,d)  eye(d,d); A [-2*sqrt(-A(1,1)) 0;0 -2*sqrt(-A(2,2))]];