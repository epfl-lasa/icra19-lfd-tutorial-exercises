function [DDX,DX,X,DX_G,Time]=simulate_modulated_system_circle(Center,Radius,X_initial,Target,A,Option)

clc
X_initial=[];
sizelim=50;
disp('Simulating the modulated dynamical system from the initial positions to the target. It might take some time, please be patient.')
x = linspace(Option.limits(1,1),Option.limits(1,2),sizelim);
lim=Option.limits;
X_d=[x      x repmat(lim(1),1,sizelim) repmat(lim(2),1,sizelim);
     repmat(lim(3),1,sizelim) repmat(lim(4),1,sizelim)    x x];
counter=1;
for i=1:size(X_d,2)
    if((Center(1)-X_d(1,i))^2+(Center(2)-X_d(2,i))^2-(Radius+Option.rho)^2>0)
        X_initial(:,counter)=X_d(:,i);
        counter=counter+1;
    end
end

for j=1:size(X_initial,2)
    j
    [DDX{j},DX{j},X{j},DX_G{j},Time{j}] = simulate(X_initial(:,j),A,Target,Center,Radius,Option);
end

end

function [DDX,DX,X,DX_G,Time]=simulate(X_initial,A,X_target,Center,Radius,option)
Deltat=option.Deltat;
delta_dx=option.delta_dx;
T=option.Tfinal;

sizeT=int64(T/Deltat);
DDX=zeros(size(X_initial,1),sizeT+1);DX=zeros(size(X_initial,1),sizeT+1);
X=zeros(size(X_initial,1),sizeT+1);
F=zeros(size(X_initial,1),sizeT+1);
DX_G=zeros(1,sizeT+1);
Time=zeros(sizeT+1,1);

counter=1;
A_2=A(3:4,3:4);
A_1=A(3:4,1:2);
X(:,1)=X_initial(1:2,1);
DX(:,1)=zeros(2,size(X_initial,2));
CONTACT=0;
rho=option.rho;

while ((counter<sizeT))
    
    Gamma=(Center(1)-X(1,counter))^2+(Center(2)-X(2,counter))^2-Radius^2;
    q1=[2*(Center(1)-X(1,counter));(Center(2)-X(2,counter)) ];
    q1=-q1/norm(q1);
    q2=[-q1(2);q1(1)];
    Q=[q1 q2;];
    Qinv=inv(Q);
        f_x=SE(DX(:,counter),X(:,counter),A_1,A_2,X_target);
    [M,~,CONTACT]=Modulation(Gamma,Radius,f_x,DX(:,counter),delta_dx,rho,q1,q2,Q,Qinv,CONTACT);
    
    if (Gamma <=0)
        F(:,counter)=(exp(-10*Gamma))*q1;
        break
    end
    DDX(:,counter+1)= M*f_x+F(:,counter);
    DX(:,counter+1)=DX(:,counter)+DDX(:,counter+1)*Deltat;
    X(:,counter+1)=X(:,counter)+DX(:,counter+1)*Deltat;
    Time(counter+1)=Time(counter)+Deltat;
    DX_G(counter+1)=q1'*DX(:,counter+1);
    if (norm(DX(:,counter+1))<0.0000001)
        break
    end
    
%     if ((rem(counter,100)==0))
%         plot(X(1,counter+1),X(2,counter+1),'.','Color',[0 0 0])
%         pause(0.0001);
%         hold on
%     end
    counter=counter+1;
end
DDX(:,counter-1:end)=[];
DX(:,counter-1:end)=[];
X(:,counter-1:end)=[];
F(:,counter-1:end)=[];
Time(counter-1:end)=[];


end
function f_x = SE(DX,X,A_1,A_2,X_target)

f_x=A_2*DX+A_1*(X-X_target);

end

function [M,Lambda,CONTACT] = Modulation(Gamma,Radius,f_x,DX,delta_dx,rho,q1,q2,Q,Qinv,Contact)
DX_G=q1'*DX;
lambda=zeros(4,1);
epsilon=5;
omega=100;
if Contact==1
    CONTACT=1;
else
    CONTACT=0;
end
rho=rho*rho+2*Radius*rho;
Deltaq_1=DX'*[2 0;0 2]*DX;
% Deltaq_1=0;
nu=-delta_dx/2;
f1=transpose(f_x)*q1/(transpose(f_x)*f_x);
f2=transpose(f_x)*q2/(transpose(f_x)*f_x);
    lambda(3) =0;
    lambda(4) =1;
if (rho<=Gamma)
  lambda(1) =((Deltaq_1-2*omega*DX_G-omega*omega*Gamma)*f1-1)*exp(epsilon*(rho-Gamma))+1;
  lambda(2) =((Deltaq_1-2*omega*DX_G-omega*omega*Gamma)*f2)*exp(epsilon*(rho-Gamma));
elseif (0<Gamma)&&(Gamma<rho)
    if (DX_G<delta_dx)
        lambda(1) =(Deltaq_1-omega*(DX_G-(delta_dx+nu)))*f1;
        lambda(2) =(Deltaq_1-omega*(DX_G-(delta_dx+nu)))*f2;
    elseif (((delta_dx<=DX_G))&&(DX_G<=0))
        lambda(1) =(Deltaq_1+(omega*omega*Gamma+nu*omega)*DX_G/delta_dx-omega*omega*Gamma)*f1;
        lambda(2) =(Deltaq_1+(omega*omega*Gamma+nu*omega)*DX_G/delta_dx-omega*omega*Gamma)*f2;
    else
        lambda(1) =(Deltaq_1-2*omega*DX_G-omega*omega*Gamma)*f1;
        lambda(2) =(Deltaq_1-2*omega*DX_G-omega*omega*Gamma)*f2;
    end
elseif (Gamma<=0)
    CONTACT=1;
    lambda(1) =(Deltaq_1-2*omega*DX_G-omega*omega*Gamma)*f1;
    lambda(2) =(Deltaq_1-2*omega*DX_G-omega*omega*Gamma)*f2;
end
Lambda=[lambda(1) lambda(2);lambda(3) lambda(4) ];
M=Q*Lambda*Qinv;

end
