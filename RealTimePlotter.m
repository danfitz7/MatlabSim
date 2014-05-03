function varargout = RealTimePlotter(varargin)
    % REALTIMEPLOTTER MATLAB code for RealTimePlotter.fig
    %      REALTIMEPLOTTER, by itself, creates a new REALTIMEPLOTTER or raises the existing
    %      singleton*.
    %
    %      H = REALTIMEPLOTTER returns the handle to a new REALTIMEPLOTTER or the handle to
    %      the existing singleton*.
    %
    %      REALTIMEPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in REALTIMEPLOTTER.M with the given input arguments.
    %
    %      REALTIMEPLOTTER('Property','Value',...) creates a new REALTIMEPLOTTER or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before RealTimePlotter_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to RealTimePlotter_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help RealTimePlotter

    % Last Modified by GUIDE v2.5 03-May-2014 13:53:52

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @RealTimePlotter_OpeningFcn, ...
                       'gui_OutputFcn',  @RealTimePlotter_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT

%SETUP
% --- Executes just before RealTimePlotter is made visible.
function RealTimePlotter_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to RealTimePlotter (see VARARGIN)

    % Choose default command line output for RealTimePlotter
    handles.output = hObject;

    % Update handles structure
    %guidata(hObject, handles);

    % UIWAIT makes RealTimePlotter wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    % START USER CODE

    %initial time variables
    global time;
    time=0;
    set(handles.periodsldr,'Value',0.01);

    % Create a timer object to fire at 1/10 sec intervals
    % Specify function handles for its start and run callbacks
    handles.timer = timer(...
        'ExecutionMode', 'fixedRate', ...       % Run timer repeatedly
        'Period', get(handles.periodsldr,'Value'), ...                        % Initial period is 1 sec.
        'TimerFcn', {@update_display,hObject}); % Specify callback function
    % Initialize slider and its readout text field
    set(handles.periodsldr,'Min',0.01,'Max',2)
    set(handles.periodsldr,'Value',get(handles.timer,'Period'))
    set(handles.slidervalue,'String',num2str(get(handles.periodsldr,'Value')))

    %%%%%%%%%%%%%%%%%%%%% CONTAINER %%%%%%%%%%%%%%%%%%%%%
    handles.containerSize=100.0;
    box=unique([0 0 0; perms([1 1 0]); perms([1 0 0]); 1 1 1], 'rows');
    box=box.*handles.containerSize;
    handles.container=box; %[0,handles.containerSize,0,handles.containerSize,0,handles.containerSize,0],Filled:FALSE, LineColor:Red);

    %plot the container on the large axis
    subplot(handles.simDisplay);                %plot on the correct axis
    set(gcf, 'CurrentAxes', handles.simDisplay);%set the axis (does this even do anything?)
    patch('Vertices',handles.container,'Faces',[1,2,4,3; 5,6,8,7; 1,2,6,5; 3,4,8,7; 1,3,7,5; 2,4,8,6],'FaceVertexCData',hsv(6),'FaceColor','flat','FaceAlpha',0.5);
    hold on

    %%%%%%%%%%%%%%%%%%%%% EFFUSION HOLE %%%%%%%%%%%%%%%%%%%%%
    holeRadius=10;
    handles.holeRadiusSquared=holeRadius^2;
    handles.holePosition=[50,50];   %on the XZ plane
    holeHeight=10;
    nHoleVerts=100;
    [x,y,z]=cylinder([holeRadius], nHoleVerts);

    x=x+handles.holePosition(1);    %translate x
    y=y+handles.holePosition(2);    %translate y
    z=z*holeHeight-holeHeight;      %scale height (z)

    %flip the hole to the X axis
    holex=z;
    holey=y;
    holez=x;  

    surf(holex,holey,holez);
    hold on

    %%%%%%%%%%%%%%%%%%%%% PARTICLES %%%%%%%%%%%%%%%%%%%%%
    handles.nParticles=1000;
    handles.temperature=1;
    handles.particleRadius=1.5;
    handles.collisionDistance=(2*handles.particleRadius)^2; %square of the sum of the radii of two particles (distance between them when they are colliding)
    handles.particleMass=1.0;
    handles.particles = handles.particleRadius+(handles.containerSize-2*handles.particleRadius)*rand(handles.nParticles,3);
    handles.particleVelocities=2*handles.temperature*rand(handles.nParticles,3)-handles.temperature;
    handles.combinations=nchoosek([1:1:handles.nParticles],2);  %1000 choose 2 oh boy!  
    
    %render as spheres
    % [SphereX,SphereY,SphereZ]=sphere(10); %make a basic n x n sphere
    % SphereX=handles.particleRadius*SphereX;
    % SphereY=handles.particleRadius*SphereY;
    % SphereZ=handles.particleRadius*SphereZ;
    % handles.sphere=[SphereX,SphereY,SphereZ];

    %handles.particleSpheresPlotHandle = surf(SphereX+handles.particles(:,1), SphereY+handle.particles(:,2), SphereZ+handles.particles(:,3));
    handles.scatterPlotHandle = scatter3(handles.simDisplay,handles.particles(:,1),handles.particles(:,2),handles.particles(:,3), 'fill');     %plot particles as 3D scatter in the main 3D axis
    %set(handles.scatterPlotHandle,'XDataSource',handles.particles(1,:),'YDataSource',handles.particles(2,:),'ZDataSource',handles.particles(3,:)); %doesn't seem to work
    
    %SIMULATION CONTROL VARIABLES
    handles.effusionEnabled=0;   %allows particles to escape through the hole when on
    handles.collisionEnabled=0;  %enables collisions between particles
    
    %%%%%%%%%%%%%%%%%%%%% initial plotting/viewport settings %%%%%%%%%%%%%%%%%%%%%

    %set axis limits
    set(handles.simDisplay,'XLim',[-50,150])
    set(handles.simDisplay,'YLim',[-50,150]);
    set(handles.simDisplay,'ZLim',[-50,150]);

    %set axis color
    set(handles.simDisplay,'XColor',[1,0,0]);
    set(handles.simDisplay,'YColor',[0,1,0]);
    set(handles.simDisplay,'ZColor',[0,0,1]);

    %kill auto-scaling
    set(handles.simDisplay,'XLimMode','manual');
    set(handles.simDisplay,'YLimMode','manual');
    set(handles.simDisplay,'ZLimMode','manual');
    set(handles.simDisplay,'DataAspectRatioMode','manual');
    set(handles.simDisplay,'DataAspectRatio',[1,1,1]);
    set(handles.simDisplay,'PlotBoxAspectRatioMode','manual');
    set(handles.simDisplay,'PlotBoxAspectRatio',[1,1,1]);

    %initial camera position
    set(handles.simDisplay,'CameraPosition',[130,125,120]);
    set(handles.simDisplay,'CameraPositionMode','manual')
    set(handles.simDisplay,'CameraTarget',[50,50,50]);
    set(handles.simDisplay,'CameraTargetMode','manual');
    set(handles.simDisplay,'Projection','perspective');
    set(handles.simDisplay,'CameraViewAngle',30);

    %PRESSURE PLOT
    %set Axis for Pressure plot
    subplot(handles.tempPlotAxis);                %plot on the correct axis
    set(gcf, 'CurrentAxes', handles.tempPlotAxis);%set the axis (does this even do anything?)
    handles.pressures=[];
    handles.times=[];
    handles.pressurePlotHandle=plot([0],[0]);%handles.tempPlotAxis,handles.times,handles.pressures);

    %N PARTICLES PLOT
    subplot(handles.nParticlesPlotAxis);                %plot on the correct axis
    set(gcf, 'CurrentAxes', handles.nParticlesPlotAxis);%set the axis (does this even do anything?)
    handles.nParticlesData=[];
    handles.nParticlePlotHandle=plot([0],[0]);
    handles.nParticlesPlotAxis.YLim.x=0;
    handles.nParticlesPlotAxis.YLim.y=handles.nParticles+1;
    
    % Update handles structure
    guidata(hObject,handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = RealTimePlotter_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

%UPDATE LOOP
%timer update
function update_display(hObject,eventdata,hfigure)
    % Timer timer1 callback, called each time timer iterates.
    
    handles = guidata(hfigure);
    axes=[1,2,3];
    
    %UPDATE TIMES
    global time;
    time=time+1;
    handles.times=[handles.times,time];
        
    %%UPDATE PARTICEL COLLISIONS
    %add particle velocities to the positions each time step
    handles.particles=handles.particles + handles.particleVelocities;
    set(handles.scatterPlotHandle,'XData',handles.particles(:,1),'YData',handles.particles(:,2),'ZData',handles.particles(:,3));
    %set(handles.particleSpheresPlotHandle, 'XData',handles.sphere(:,1)+handles.particles(:,1), 'YData',handles.sphere(:,2)+handles.particles(:,2), 'ZData',handles.sphere(:,2)+handles.particles(:,2));

    %%ESCAPE BY EFFUSION THROUGH HOLE
    if (handles.effusionEnabled==1)
        effusingParticlesIndeces=find((handles.particles(:,1)-handles.particleRadius<0) & (((handles.particles(:,2)-handles.holePosition(1)).^2+(handles.particles(:,3)-handles.holePosition(2)).^2)<handles.holeRadiusSquared));
        handles.particles=removerows(handles.particles,effusingParticlesIndeces);
        handles.particleVelocities=removerows(handles.particleVelocities,effusingParticlesIndeces);
        handles.nParticles=size(handles.particles,1);
    end
    
    %Update nParticles graph
    handles.nParticlesData=[handles.nParticlesData,handles.nParticles];
    set(handles.nParticlePlotHandle,'XData',handles.times,'Ydata',handles.nParticlesData);
    
    %%COLLISION WITH CONTAINER
    %colliding=@(a,i)  (handles.particles(a,i)>containerSize);
    wallPressure=0;
    for axis=axes
        collidingParticleIndeces=ind2sub(size(handles.particles(:,axis)),find((handles.particles(:,axis)+handles.particleRadius>handles.containerSize) | (handles.particles(:,axis)-handles.particleRadius<0)));
        %wallPressure=wallPressure+sum(abs(handles.particleVelocities(collidingParticleIndeces,axis)));
        handles.particleVelocities(collidingParticleIndeces,axis) = -1.*handles.particleVelocities(collidingParticleIndeces,axis);   %+X
        handles.combinations=nchoosek([1:1:handles.nParticles],2);  %1000 choose 2 oh boy!  
    end

    %INTERPARTICLE COLLISIONS
    particleCombinationSquareDistances=sum(handles.particles((handles.combinations(:,1))-handles.particles(handles.combinations(:,2))).^2,2); %for each combination of particles, find the sum of the squares of the differences of the elements of their position vectors (the distance between them squared)
    collidingPairs=find(particleCombinationSquareDistances<=handles.collisionDistance);
    collidingPairIndeces=handles.combinations(ind2sub(size(handles.combinations), collidingPairs),:);   %a list of pairs of indeces for the particles and particleVelcotity arrays of collising particles
    
    collidingIndeces=unique(collidingPairIndeces(:));
    handles.particles=removerows(handles.particles,collidingIndeces);
    handles.particleVelocities=removerows(handles.particleVelocities,collidingIndeces);
    handles.nParticles=size(handles.particles,1);
    
    %UPDATE PRESSURES
    handles.pressureConstant = handles.particleMass/handles.containerSize; %F_AVG = (m N v_avg)/L
    curPressure=handles.pressureConstant*sumsqr(handles.particleVelocities);  %mean(handles.particleVelocities(:,1).^2+handles.particleVelocities(:,2).^2+handles.particleVelocities(:,3).^2)
    handles.pressures=[handles.pressures,curPressure];%[handles.pressures,wallPressure];
    set(handles.pressurePlotHandle,'XData',handles.times,'Ydata',handles.pressures);

    % Update handles structure
    guidata(hfigure,handles);

    %refreshdata(handles.simDisplay);   %unneccisary?
    drawnow;
    % END USER CODE
end

% --- Executes on button press in startbtn.
function startbtn_Callback(hObject, eventdata, handles)
    % hObject    handle to startbtn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % START USER CODE
    % Only start timer if it is not running
    if strcmp(get(handles.timer, 'Running'), 'off');
        start(handles.timer);
        disp('STARTING...');
    end
end

% --- Executes on button press in stopbtn.
function stopbtn_Callback(hObject, eventdata, handles)
    % hObject    handle to stopbtn (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % START USER CODE
    % Only stop timer if it is running
    if strcmp(get(handles.timer, 'Running'), 'on')
        stop(handles.timer);
        disp('Paused.');
    end
end

% --- Executes on slider movement.
function periodsldr_Callback(hObject, eventdata, handles)
    % hObject    handle to periodsldr (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'Value') returns position of slider
    %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    % START USER CODE
    % Read the slider value
    period = get(handles.periodsldr,'Value');

    % Timers need the precision of periods to be greater than about
    % 1 millisecond, so truncate the value returned by the slider
    period = period - mod(period,.01);

    % Set slider readout to show its value
    set(handles.slidervalue,'String',num2str(period))

    % If timer is on, stop it, reset the period, and start it again.
    if strcmp(get(handles.timer, 'Running'), 'on')
        stop(handles.timer);
        set(handles.timer,'Period',period)
        start(handles.timer)
    else               % If timer is stopped, reset its period only.
        set(handles.timer,'Period',period)
    end

    disp('SLIDER MOVED');
    % END USER CODE
end

% --- Executes on button press in EffusionChechBox.
function EffusionChechBox_Callback(hObject, eventdata, handles)
    % hObject    handle to EffusionChechBox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of EffusionChechBox
    handles.effusionEnabled = get(hObject,'Value');
    guidata(hObject,handles);
end

% --- Executes on button press in ParticleCollisionCheckbox.
function ParticleCollisionCheckbox_Callback(hObject, eventdata, handles)
    % hObject    handle to ParticleCollisionCheckbox (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hint: get(hObject,'Value') returns toggle state of ParticleCollisionCheckbox
    handles.collisionEnabled = get(hObject,'Value');
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function periodsldr_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to periodsldr (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
end

% --- Executes during object creation, after setting all properties.
function simDisplay_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to simDisplay (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: place code in OpeningFcn to populate simDisplay
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
    % hObject    handle to figure1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % START USER CODE
    % Necessary to provide this function to prevent timer callback
    % from causing an error after GUI code stops executing.
    % Before exiting, if the timer is running, stop it.
    if strcmp(get(handles.timer, 'Running'), 'on')
        stop(handles.timer);
    end
    % Destroy timer
    delete(handles.timer)
    % END USER CODE

    % Hint: delete(hObject) closes the figure
    delete(hObject);
end
