function [DDX,DX,X,F,Time]=simulate_modulated_system(A,N_X,Poly,X_initial,X_target,X_free,X_C,X_L,option)



disp('Simulating the modulated dynamical system from the initial positions to the target. It might take some time, please be patient.')

for j=1:size(X_initial,2)
    [DDX{j},DX{j},X{j},F{j},Time{j}] = simulate(X_initial(:,j),A,X_target,X_free,X_C,X_L,Poly,N_X,option);
end

end

function [DDX,DX,X,F,Time]=simulate(X_initial,A,X_target,X_free,X_C,X_L,Poly,N_X,option)
Deltat=option.Deltat;
% F_d=option.F_d;
delta_dx=option.delta_dx;
T=option.Tfinal;

Handle_sign=sign(-X_free(2,1)+Poly(1)*X_free(1,1)+Poly(2));
sizeT=int64(T/Deltat);
DDX=zeros(size(X_initial,1)/2,sizeT+1);DX=zeros(size(X_initial,1)/2,sizeT+1);
X=zeros(size(X_initial,1)/2,sizeT+1);
F=zeros(size(X_initial,1)/2,sizeT+1);
Time=zeros(sizeT+1,1);

q2=[-N_X(2);N_X(1)];
Q=[N_X q2];
Qinv=inv(Q);

counter=1;
A_2=A(3:4,3:4);
A_1=A(3:4,1:2);



X(:,1)=X_initial(1:2,1);
DX(:,1)=10*X_initial(3:4,1);
% X_mu=(X_C+X_L)/2;
% signal_i=eye(2);
CONTACT=0;
rho=option.rho;
% Limits=option.limits;
% screensize = get( 0, 'Screensize' );
% fig = figure();
% subplot1 = subplot(1,1,1);
% set(fig,'Position',screensize)
% plot_Wall_counters(subplot1,N_X,X_free,X_C,X_L,rho,Limits)
% h2=plot(X_target(1,1),X_target(2,1),'MarkerFaceColor',[0 0 1],...
%     'MarkerEdgeColor','none',...
%     'MarkerSize',30,...
%     'Marker','pentagram',...
%     'LineWidth',5,...
%     'LineStyle','none');
% hold on
% h1=plot(X_initial(1,:),X_initial(2,:),...
%     'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
%     'MarkerEdgeColor','none',...
%     'MarkerSize',30,...
%     'Marker','hexagram',...
%     'LineWidth',5,...
%     'LineStyle','none');
% X_wall(1,:) = linspace(Limits(1,1),Limits(1,2),10^5);
% X_wall(2,:) = Poly(1)*X_wall(1,:)+Poly(2);
% h3=plot(X_wall(1,:),X_wall(2,:),'LineWidth',4,...
%     'LineStyle','--',...
%     'Color',[0 0 0]);
% hold on
% h4=plot(X_C(1,1),X_C(2,1),...
%     'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
%     'MarkerSize',30,...
%     'Marker','^',...
%     'LineStyle','none');
% h5=plot(X_L(1,1),X_L(2,1),...
%     'MarkerFaceColor',[1 0 0],...
%     'MarkerSize',30,...
%     'Marker','v',...
%     'LineStyle','none');
Wall_Base=N_X'*X_C;
Handle_sign=sign(N_X'*X_free-Wall_Base);
while ((counter<sizeT))
    
    Q_2=(X_L-X_C);
    Gamma=Handle_sign*(N_X'*X(:,counter)-Wall_Base)+(rho-(X_L-X_C)'*(X_L-X(:,counter)))*exp(-option.kamma_slider*(X_L-X(:,counter))'*(X_L-X(:,counter)));
    %     Gamma=((5-exp(-transpose(X(:,counter)-X_mu)*signal_i*(X(:,counter)-X_mu)))/(1+exp(-transpose(X(:,counter)-X_mu)*signal_i*(X(:,counter)-X_mu)))-2)...
    %         *(0.1*Handle_sign*(-X(2,counter)+Poly(1)*X(1,counter)+Poly(2))+exp(-10000*transpose(X_L-X(:,counter))*Q_2));
    %    Gamma=((5-exp(-transpose(X(:,counter)-X_mu)*signal_i*(X(:,counter)-X_mu)))/(1+exp(-transpose(X(:,counter)-X_mu)*signal_i*(X(:,counter)-X_mu)))-2)...
    %         *(0.1*+exp(-10000*transpose(X_L-X(:,counter))*Q_2));
    
    if CONTACT
        f_x=SE(DX(:,counter),X(:,counter),100*A_1,10*A_2,X_target);
    else
        f_x=SE(DX(:,counter),X(:,counter),A_1,A_2,X_target);
    end
    [M,~,CONTACT]=Modulation_2(Gamma,Wall_Base,f_x,DX(:,counter),X(:,counter),X_C,X_L,delta_dx,rho,N_X,q2,Q,Qinv,CONTACT);
    
    if (Handle_sign*(N_X'*X(:,counter)-Wall_Base) <=0)
        F(:,counter)=(exp(-(Handle_sign*(-X(2,counter)+Poly(1)*X(1,counter)+Poly(2))))-1)*N_X;
        X(:,counter)=(X(:,counter)-X_L)'*(Q_2)*Q_2/(norm(Q_2)^2)+X_L;
    end
    DDX(:,counter+1)= M*f_x+F(:,counter);
    DX(:,counter+1)=DX(:,counter)+DDX(:,counter+1)*Deltat;
    X(:,counter+1)=X(:,counter)+DX(:,counter+1)*Deltat;
    Time(counter+1)=Time(counter)+Deltat;
    if (norm(X(:,counter+1)-X_target)<0.1)
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

f_x=2*A_2*DX+A_1*(X-X_target);

end

function [M,Lambda,CONTACT] = Modulation(Gamma,f_x,DX,X,X_C,F_d,delta_dx,rho,N_X,q2,Q,Qinv,Contact)
epsilon=0.00000001;
DX_G=N_X'*DX;
Lambda=zeros(2,2);
if Contact==1
    CONTACT=1;
else
    CONTACT=0;
end
if (1<=Gamma)
    Lambda(1,1) = (1-exp(-(1/epsilon)*(Gamma-1)))/(1+exp(-(1/epsilon)*(Gamma-1)));
end
if (((DX_G<delta_dx)&&((0<Gamma)&&(Gamma<1)))&&( Contact==0))
    Lambda(1,1) = 10*(1-Gamma)*(delta_dx-DX_G)/(abs(Gamma)*N_X'*f_x);
end
if (((DX_G<0)&&((delta_dx<=DX_G)))&&((0<Gamma)&&(Gamma<1)))
    Lambda(1,1) =-(1-Gamma)*F_d*exp(-Gamma/epsilon);
end
if ((0<=DX_G)&&((0<Gamma)&&(Gamma<1)))
    Lambda(1,1) =  -(1-Gamma)*(DX_G+exp(-Gamma/epsilon))*F_d/(N_X'*f_x);
end
if ((Gamma<=0))
    Lambda(1,1) =  -F_d/(N_X'*f_x);
end
if (1<=Gamma)
    Lambda(2,2)=Lambda(1,1);
elseif ((0<Gamma)&&(Gamma<1)&&( Contact==0))
    Lambda(2,2)= 10*(1-Gamma)*(-(1/(Gamma+0.001))*(q2'*DX*N_X'*(X-X_C)-N_X'*DX*q2'*(X-X_C)))/(N_X'*(X-X_C)*q2'*f_x);
elseif ((Gamma<=0)&&( Contact==0))
    CONTACT=1;
    Lambda(2,2)=1;
elseif ( Contact==1)
    Lambda(2,2)=1;
end
%     if norm(Lambda)>1000
%         keyboard
%     end
M=Q*Lambda*Qinv;

end

function [M,Lambda,CONTACT] = Modulation_2(Gamma,Wall_Base,f_x,DX,X,X_c,X_l,delta_dx,rho,N_X,q2,Q,Qinv,Contact)
DX_G=N_X'*DX;
lambda=zeros(4,1);
epsilon=50;
if Contact==1
    CONTACT=1;
    omega=1;
else
    CONTACT=0;
    omega=0.01;
end

nu=0.05;
f1=transpose(f_x)*N_X/(transpose(f_x)*f_x);
f2=transpose(f_x)*q2/(transpose(f_x)*f_x);
if (rho<=Gamma)
    lambda(1) =((-2*omega^(-1)*transpose(N_X)*DX-omega^(-2)*transpose(N_X)*(X-X_c))*f1-1)*exp(epsilon*(rho-Gamma))+1;
    lambda(2) =((-2*omega^(-1)*transpose(N_X)*DX-omega^(-2)*transpose(N_X)*(X-X_c))*f2)*exp(epsilon*(rho-Gamma));
    lambda(3) =((-2*omega^(-1)*transpose(q2)*DX-omega^(-2)*transpose(q2)*(X-X_c))*f1)*exp(epsilon*(rho-Gamma));
    lambda(4) =((-2*omega^(-1)*transpose(q2)*DX-omega^(-2)*transpose(q2)*(X-X_c))*f2-1)*exp(epsilon*(rho-Gamma))+1;
elseif (0<transpose(N_X)*X-Wall_Base)&&(Gamma<rho)
    if (transpose(N_X)*DX<delta_dx)
        lambda(1) =-omega^(-1)*(transpose(N_X)*DX-(delta_dx+nu))*f1;
        lambda(2) =-omega^(-1)*(transpose(N_X)*DX-(delta_dx+nu))*f2;
    elseif (((delta_dx<=transpose(N_X)*DX))&&(transpose(N_X)*DX<=0))
        lambda(1) =0;
        lambda(2) =0;
    else
        lambda(1) =((-2*omega^(-1)*transpose(N_X)*DX-omega^(-2)*transpose(N_X)*(X-X_c))*f1);
        lambda(2) =((-2*omega^(-1)*transpose(N_X)*DX-omega^(-2)*transpose(N_X)*(X-X_c))*f2);
    end
    lambda(3) =((-2*omega^(-1)*transpose(q2)*DX-omega^(-2)*transpose(q2)*(X-X_c))*f1);
    lambda(4) =((-2*omega^(-1)*transpose(q2)*DX-omega^(-2)*transpose(q2)*(X-X_c))*f2);
elseif (transpose(N_X)*X-Wall_Base<=0)
    CONTACT=1;
    omega=1;
    lambda(1) =((-2*omega^(-1)*transpose(N_X)*DX-omega^(-2)*transpose(N_X)*(X-(2*X_l-X_c)))*f1);
    lambda(2) =((-2*omega^(-1)*transpose(N_X)*DX-omega^(-2)*transpose(N_X)*(X-(2*X_l-X_c)))*f2);
    lambda(3) =((-2*omega^(-1)*transpose(q2)*DX-omega^(-2)*transpose(q2)*(X-(2*X_l-X_c)))*f1);
    lambda(4) =((-2*omega^(-1)*transpose(q2)*DX-omega^(-2)*transpose(q2)*(X-(2*X_l-X_c)))*f2);
end
Lambda=[lambda(1) lambda(2);lambda(3) lambda(4) ];
M=Q*Lambda*Qinv;
end
