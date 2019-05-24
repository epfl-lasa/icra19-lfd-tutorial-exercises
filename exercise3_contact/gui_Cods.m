function varargout = gui_Cods(varargin)
% GUI_CODS MATLAB code for gui_Cods.fig
%      GUI_CODS, by itself, creates a new GUI_CODS or raises the existing
%      singleton*.
%
%      H = GUI_CODS returns the handle to a new GUI_CODS or the handle to
%      the existing singleton*.
%
%      GUI_CODS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CODS.M with the given input arguments.
%
%      GUI_CODS('Property','Value',...) creates a new GUI_CODS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_Cods_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_Cods_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_Cods

% Last Modified by GUIDE v2.5 28-May-2018 10:57:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_Cods_OpeningFcn, ...
    'gui_OutputFcn',  @gui_Cods_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_Cods is made visible.
function gui_Cods_OpeningFcn(hObject, eventdata, handles, varargin)
addpath(genpath('CoDS'))
setup_CoDs;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_Cods (see VARARGIN)
cla
legend('off');
handles.uipanel2.Visible='off';
handles.simulate.Visible='off';
% Choose default command line output for gui_Cods
handles.output = hObject;
handles.Deltat=0.001;
handles.check=0;
handles.check_rho=0;
handles.delta_dx=-0.2;
handles.Tfinal=10;
handles.animation=0;
handles.Onsurface=0;
handles.kamma_slider=1;
handles.rho=get(handles.Rho_slider,'Value');
handles.axes1.XLim=[-5 5];
handles.axes1.YLim=[-5 5];
% Update handles structure
% screensize = get( 0, 'Screensize' );
%  set(handles.axes1,'Position',screensize)
guidata(hObject, handles);

% UIWAIT makes gui_Cods wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.Rho_value,'String',get(handles.Rho_slider,'Value'));
grid on
box on



% --- Outputs from this function are returned to the command line.
function varargout = gui_Cods_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function Rho_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Rho_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Rho_value,'String',get(handles.Rho_slider,'Value'));
handles.rho=get(handles.Rho_slider,'Value');
guidata(hObject, handles);
cla
grid on
box on

X1=[ handles.axes1.XLim(1):0.05: handles.axes1.XLim(2)];
X2=[ handles.axes1.YLim(1):0.05: handles.axes1.YLim(2)]';

X1=repmat(X1,size(X2,1),1);
X2=repmat(X2,1,size(X1,2));
Walla=zeros(size(X2));

Wall_Base=handles.N_x'*handles.X_C;
Handle_sign=sign(handles.N_x'*handles.X_free-Wall_Base);
for ii=1:size(X2,1),
    for jj=1:size(X2,2)
        XX=[X1(ii,jj);X2(ii,jj)];
        Walla(ii,jj)=Handle_sign*(handles.N_x'*XX-Wall_Base)+...
            (handles.rho-(handles.X_L-handles.X_C)'*(handles.X_L-XX))*exp(-handles.kamma_slider*(handles.X_L-XX)'*(handles.X_L-XX));
        if  handles.rho<(Walla(ii,jj))
            Walla(ii,jj)=handles.rho;
        end
        if Walla(ii,jj)<-2
            Walla(ii,jj)=-2;
        end
    end
end
clim=[-2 handles.rho];
hold on
contourf(X1,X2,Walla);
colormap(hot);
hold on
X = linspace( handles.limits(1), handles.limits(2),100);
Y = handles.Poly(1)*X+handles.Poly(2);
h3=plot(X,Y,'LineWidth',4,'LineStyle','--','Color',[0 0 0]);
hold on
h2=plot(handles.X_target(1,1),handles.X_target(2,1),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',30,'Marker','pentagram','LineWidth',5,'LineStyle','none');
h1=plot(handles.X_initial(1,:),handles.X_initial(2,:),'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],'MarkerEdgeColor','none','MarkerSize',30,'Marker','hexagram','LineWidth',5,'LineStyle','none');
h4=plot(handles.X_C(1,1),handles.X_C(2,1),'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],'MarkerSize',30,'Marker','^','LineStyle','none');
if (handles.Onsurface==0)
    h5=plot(handles.X_L(1,1),handles.X_L(2,1),...
        'MarkerFaceColor',[1 0 0],...
        'MarkerSize',30,...
        'Marker','v',...
        'LineStyle','none');
end
% ylim([app.Option.limits(3) app.Option.limits(4)]);
% xlim([app.Option.limits(1) app.Option.limits(2)]);
grid(handles.axes1,'on')
grid on
box on


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Rho_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rho_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function X_axis_Callback(hObject, eventdata, handles)
% hObject    handle to X_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.axes1.XLim=[-get(handles.X_axis,'Value') get(handles.X_axis,'Value')];

% --- Executes during object creation, after setting all properties.
function X_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to X_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Y_axis_Callback(hObject, eventdata, handles)
% hObject    handle to Y_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.axes1.YLim=[-get(handles.Y_axis,'Value') get(handles.Y_axis,'Value')];

% --- Executes during object creation, after setting all properties.
function Y_axis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Y_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Animation.
function Animation_Callback(hObject, eventdata, handles)
% hObject    handle to Animation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.animation=get(hObject,'Value');
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Animation



function velocity_Callback(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.delta_dx=str2double(get(hObject,'String'));
if handles.delta_dx<-10
    handles.delta_dx=-10;
    set(hObject,'String','-10');
end
if handles.delta_dx>-0.001
    handles.delta_dx=-0.001;
    set(hObject,'String','-0.001');
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of velocity as text
%        str2double(get(hObject,'String')) returns contents of velocity as a double


% --- Executes during object creation, after setting all properties.
function velocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_time_Callback(hObject, eventdata, handles)
% hObject    handle to step_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Deltat=str2double(get(hObject,'String'));
if handles.Deltat<0.00001
    handles.Deltat=0.00001;
    set(hObject,'String','0.00001');
end
if handles.Deltat>0.01
    handles.Deltat=0.01;
    set(hObject,'String','0.01');
end
guidata(hObject, handles);
% handles
% Hints: get(hObject,'String') returns contents of step_time as text
%        str2double(get(hObject,'String')) returns contents of step_time as a double


% --- Executes during object creation, after setting all properties.
function step_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_Callback(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Tfinal=str2double(get(hObject,'String'));
if handles.Tfinal<1
    handles.Tfinal=1;
    set(hObject,'String','1');
end
if handles.Tfinal>100
    handles.Tfinal=100;
    set(hObject,'String','100');
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of time as text
%        str2double(get(hObject,'String')) returns contents of time as a double


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on button press in Off_surface.
function Off_surface_Callback(hObject, eventdata, handles)
% hObject    handle to Off_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uipanel2.Visible='off';
handles.simulate.Visible='off';
handles.N_x=[];
handles.Poly=[];
handles.X_target=[];
handles.check=[];
handles.A=[];
handles.X_initial=[];
handles.X_C=[];
handles.X_L=[];

cla
legend('off');
handles.limits= [ handles.axes1.XLim handles.axes1.YLim];
handles.F_d=-5;
handles.Onsurface=0;
disp('Draw the contact surface')
[handles.N_x,handles.Poly,handles.X_target,handles.check]=Construct_the_surface(handles);
if handles.check==0
    error('Program exit')
end
disp('Draw some motions, make sure that it ends up at the target point and it goes through the contact surface !')
[handles.A,handles.X_initial,handles.check]=Construct_the_dynamcial_system(handles.Poly,handles.X_target,handles.X_target,handles);
if handles.check==0
    error('Program exit')
end
[handles.X_C,handles.X_L]=Select_the_contact_point(handles.Poly,handles.X_target,handles.X_initial,handles);
handles.X_free=handles.X_target;
handles.Play_with_rho=1;
handles.uipanel2.Visible='on';
handles.simulate.Visible='on';

guidata(hObject, handles);








% --- Executes on button press in On_surface.
function On_surface_Callback(hObject, eventdata, handles)
% hObject    handle to On_surface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uipanel2.Visible='off';
handles.simulate.Visible='off';
handles.N_x=[];
handles.Poly=[];
handles.X_target=[];
handles.check=[];
handles.A=[];
handles.X_initial=[];
handles.X_C=[];
handles.X_L=[];
cla
handles.limits= [ handles.axes1.XLim handles.axes1.YLim];
handles.F_d=-5;
handles.Onsurface=1;
disp('Draw the contact surface')
[handles.N_x,handles.Poly,handles.X_target,handles.X_free,handles.check]=Construct_the_surface_for_on_surface(handles);
if handles.check==0
    error('Program exit')
end
disp('Draw some motions, make sure that it ends up at the target point and it goes through the contact surface !')
[handles.A,handles.X_initial,handles.Option.check]=Construct_the_dynamcial_system(handles.Poly,handles.X_target,handles.X_free,handles);
if handles.check==0
    error('Program exit')
end
[handles.X_C,handles.X_L]=Select_the_contact_point_on_surface(handles.Poly,handles.X_target,handles.X_initial,handles);



handles.uipanel2.Visible='on';
handles.simulate.Visible='on';
guidata(hObject, handles);

% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[DDX_modulated,DX_modulated,X_modulated,F_modulated,Time_modulated]=  simulate_modulated_system(handles.A,handles.N_x,handles.Poly,handles.X_initial,handles.X_target,handles.X_free,handles.X_C,handles.X_L,handles);
plot_the_simualtions(DDX_modulated,DX_modulated,X_modulated,F_modulated,Time_modulated,handles.Poly,handles.X_initial,handles.X_target,handles.X_free,handles.X_C,handles.X_L,handles.N_x,handles);


% --- Executes on button press in stop_record.
function stop_record_Callback(hObject, eventdata, handles)
% hObject    handle to stop_record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.stop_record.Value=1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on slider movement.
function kamma_slider_Callback(hObject, eventdata, handles)
% hObject    handle to kamma_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% set(handles.kamma_value,'String',2);
handles.kamma_slider=get(hObject,'Value');
set(handles.kamma_value,'String',get(hObject,'Value'));
handles.kamma=get(handles.kamma_value,'Value');
guidata(hObject, handles);
cla
grid on
box on

X1=[ handles.axes1.XLim(1):0.05: handles.axes1.XLim(2)];
X2=[ handles.axes1.YLim(1):0.05: handles.axes1.YLim(2)]';

X1=repmat(X1,size(X2,1),1);
X2=repmat(X2,1,size(X1,2));
Walla=zeros(size(X2));

Wall_Base=handles.N_x'*handles.X_C;
Handle_sign=sign(handles.N_x'*handles.X_free-Wall_Base);
for ii=1:size(X2,1),
    for jj=1:size(X2,2)
        XX=[X1(ii,jj);X2(ii,jj)];
        Walla(ii,jj)=Handle_sign*(handles.N_x'*XX-Wall_Base)+...
            (handles.rho-(handles.X_L-handles.X_C)'*(handles.X_L-XX))*exp(-handles.kamma_slider*(handles.X_L-XX)'*(handles.X_L-XX));
        if  handles.rho<(Walla(ii,jj))
            Walla(ii,jj)=handles.rho;
        end
        if Walla(ii,jj)<-2
            Walla(ii,jj)=-2;
        end
    end
end
clim=[-2 handles.rho];
hold on
contourf(X1,X2,Walla);
colormap(hot);
hold on
X = linspace( handles.limits(1), handles.limits(2),100);
Y = handles.Poly(1)*X+handles.Poly(2);
h3=plot(X,Y,'LineWidth',4,'LineStyle','--','Color',[0 0 0]);
hold on
h2=plot(handles.X_target(1,1),handles.X_target(2,1),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor','none','MarkerSize',30,'Marker','pentagram','LineWidth',5,'LineStyle','none');
h1=plot(handles.X_initial(1,:),handles.X_initial(2,:),'MarkerFaceColor',[0.466666668653488 0.674509823322296 0.18823529779911],'MarkerEdgeColor','none','MarkerSize',30,'Marker','hexagram','LineWidth',5,'LineStyle','none');
h4=plot(handles.X_C(1,1),handles.X_C(2,1),'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],'MarkerSize',30,'Marker','^','LineStyle','none');
if (handles.Onsurface==0)
    h5=plot(handles.X_L(1,1),handles.X_L(2,1),...
        'MarkerFaceColor',[1 0 0],...
        'MarkerSize',30,...
        'Marker','v',...
        'LineStyle','none');
end

% --- Executes during object creation, after setting all properties.
function kamma_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kamma_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
