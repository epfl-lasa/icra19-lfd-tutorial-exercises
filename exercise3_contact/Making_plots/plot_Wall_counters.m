function plot_Wall_counters(subplot1,N_X,X_free,X_C,X_L,rho,Limits)

hold(subplot1,'on');
box(subplot1,'on');
grid(subplot1,'on');
set(subplot1,'FontSize',20);
xlabel('X(1)','Interpreter','latex');
ylabel('X(2)','Interpreter','latex');
ylim([Limits(3) Limits(4)]);
xlim([Limits(1) Limits(2)]);

X1=[Limits(1):0.05:Limits(2)];
X2=[Limits(3):0.05:Limits(4)]';
X1=repmat(X1,size(X2,1),1);
X2=repmat(X2,1,size(X1,2));


Walla=zeros(size(X2));
Wall_Base=N_X'*X_C;

Handle_sign=sign(N_X'*X_free-Wall_Base);
for ii=1:size(X2,1),
    for jj=1:size(X2,2)
        XX=[X1(ii,jj);X2(ii,jj)];
        Walla(ii,jj)=Handle_sign*(N_X'*XX-Wall_Base)+...
            (rho-(X_L-X_C)'*(X_L-XX))*exp(-2*(X_L-XX)'*(X_L-XX));
        if  rho<(Walla(ii,jj))
            Walla(ii,jj)=rho;
        end
    end
end
clim=[-2 rho];
contourf(X1,X2,Walla)
colormap(hot)
set(subplot1,'CLim',clim);
set(subplot1,'BoxStyle','full','CLim',clim,'FontSize',20,'Layer','top',...
    'TickLabelInterpreter','latex');
colorbar('peer',subplot1);
hold on
