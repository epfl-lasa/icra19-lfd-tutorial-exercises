function [data,x_obs_dem_initial] = generate_mouse_data(limits,derivative,filter,Title)
% GENERATE_MOUSE_DATA(NTH_ORDER, N_DOWNSAMPLE) request the user to give
% demonstrations of a trajectories in a 2D workspace using the mouse cursor
% The data is stored in an [x ; dx/dt] structure
% The data isdownsampled by N_DOWNSAMPLE samples

%   # Authors: Klas Kronander, Wissam Bejjani, Jose Medina and Seyed Sina
%   Mirarzavi Salehian
%   # EPFL, LASA laboratory
%   # Email: sina.mirrazavi@epfl.com

%% Drawing plot
% fig = figure();
% view([0 90]);

disp('Draw some trajectories with the mouse on the figure.')
% axes('Parent',fig);
title(Title);
 ylim([limits(3) limits(4)]);
 xlim([limits(1) limits(2)]);

box on;
grid  on;
hold on;
% ylabel('X(2)','Interpreter','latex');
% xlabel('X(1)','Interpreter','latex');
delete_trace = 0;

% to store the data
X = [];
% flag for signaling that the demonstration has ended
demonstration_index = 0;
demonstration_index_monitor = 0;

% select our figure as gcf
% figure(fig);
hold on
% disable any figure modes
zoom off
rotate3d off
pan off
brush off
datacursormode off

set(gcf,'WindowButtonDownFcn',@(h,e)button_clicked(h,e));
set(gcf,'WindowButtonUpFcn',[]);
set(gcf,'WindowButtonMotionFcn',[]);
set(gcf,'Pointer','circle');
hp = gobjects(0);
%       'position',[1350 432 300 100]
  stop_btn = uicontrol('style','pushbutton','String', 'stop recording','Callback',@stop_recording,'Units','normalized', ...
      'position',[0.789 0.398 0.16 0.093], ...
      'UserData', 1);
% option.stop_record.Visible='on';
% wait until demonstration is finished
while( (get(stop_btn, 'UserData') == 1));
    pause(0.1);
    if demonstration_index ~= demonstration_index_monitor;
        x_obs{demonstration_index} = X;
        X = [];
        demonstration_index_monitor = demonstration_index;
        set(gcf,'WindowButtonDownFcn',@(h,e)button_clicked(h,e));
        set(gcf,'WindowButtonUpFcn',[]);
        set(gcf,'WindowButtonMotionFcn',[]);
        set(gcf,'Pointer','circle');
    end
end
n_demonstrations = demonstration_index_monitor;
delete(stop_btn);
set(gcf,'Pointer','arrow');
%% Savitzky-Golay filter and derivatives
%   x :             input data size (time, dimension)
%   dt :            sample time
%   nth_order :     max order of the derivatives
%   n_polynomial :  Order of polynomial fit
%   window_size :   Window length for the filter

if filter ==1
    x_obs_dem_initial=zeros(4,n_demonstrations);
    if derivative==1
        data = [];
        for dem = 1:n_demonstrations
            x_obs_dem = x_obs{dem}(1:2,:)';
            x_obs_dem_initial(1:2,dem)=x_obs_dem(1,:)';
            dt = mean(diff(x_obs{dem}(3,:)')); % Average sample time (Third dimension contains time)
            dx_nth = sgolay_time_derivatives(x_obs_dem, dt, 2, 3, min(9,2*floor(size(x_obs_dem,1)/2)+1));
            x_obs_dem_initial(3:4,dem)=dx_nth(1,:,2)';
            data = [data [dx_nth(:,:,1),dx_nth(:,:,2),dx_nth(:,:,2),dx_nth(:,:,3)]'];
        end
        
    else
        data = [];
        for dem = 1:n_demonstrations
            x_obs_dem = x_obs{dem}(1:2,:)';
            x_obs_dem_initial(1:2,dem)=x_obs_dem(1,:)';
            dt = mean(diff(x_obs{dem}(3,:)')); % Average sample time (Third dimension contains time)
            dx_nth = sgolay_time_derivatives(x_obs_dem, dt, 2, 3, min(9,2*floor(size(x_obs_dem,1)/2)+1));
             x_obs_dem_initial(3:4,dem)=dx_nth(1,:,2)';
            data = [data [dx_nth(:,:,1)]'];
        end
    end
else
    data = [];
    x_obs_dem_initial=zeros(2,n_demonstrations);
    for dem = 1:n_demonstrations
        x_obs_dem = x_obs{dem}(1:2,:)';
        x_obs_dem_initial(:,dem)=x_obs_dem(1,:)';
        data = [data [x_obs_dem]'];
    end
end

return


%% Functions for data capture
    function stop_recording(ObjectS, ~)
        set(ObjectS, 'UserData', 0);
    end

    function ret = button_clicked(~,~)
        if(strcmp(get(gcf,'SelectionType'),'normal'));
            start_demonstration();
        end
    end

    function ret = start_demonstration()
        disp('Started demonstration');
        set(gcf,'WindowButtonUpFcn',@stop_demonstration);
        set(gcf,'WindowButtonMotionFcn',@record_current_point);
        ret = 1;
        tic;
    end

    function ret = stop_demonstration(~,~)
        disp('Stopped demonstration. Press stop recording in the figure if you have enough demonstrations.');
        set(gcf,'WindowButtonMotionFcn',[]);
        set(gcf,'WindowButtonUpFcn',[]);
        set(gcf,'WindowButtonDownFcn',[]);
        if(delete_trace);
            delete(hp);
        end
        demonstration_index = demonstration_index + 1;
    end

    function ret = record_current_point(~,~)
        x = get(gca,'Currentpoint');
        x = x(1,1:2)';
        x = [x;toc];
        X = [X, x];
        hp = [hp, plot(x(1),x(2),'r.','markersize',20)];
    end
end


