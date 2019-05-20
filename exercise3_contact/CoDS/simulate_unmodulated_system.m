function simulate_unmodulated_system(A,Dem_i)

options = odeset('RelTol',1e-4,'AbsTol',1e-4);
for j=1:size(Dem_i,2)
    Xrobot_object_initial=Dem_i(:,j);
%  Xrobot_object_initial=[Dem_i(1:2,j)];
[~,x{j}] = ode45(@(t,xi_r) SE(t,xi_r,A),[0:0.001:50],Xrobot_object_initial,options);
end

for j=1:size(Dem_i,2)
 plot(x{j}(1,1),x{j}(1,2),'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
    'MarkerEdgeColor',[0.466666668653488 0.674509823322296 0.18823529779911],...
    'MarkerSize',20,...
    'Marker','hexagram',...
    'LineStyle','none','DisplayName','Intiial Position');
plot(x{j}(:,1),x{j}(:,2),'DisplayName','Generated trajectory','LineWidth',1,'Color',[0 0 0]);
hold on


end




 
 
function Dx_main = SE(t,xi_r,A)
Dx_main=  A* xi_r;



