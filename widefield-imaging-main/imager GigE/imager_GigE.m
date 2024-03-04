function varargout = imager_GigE(varargin)
% IMAGER_GIGE MATLAB code for imager_GigE.fig
%      IMAGER_GIGE, by itself, creates a new IMAGER_GIGE or raises the existing
%      singleton*.
%
%      H = IMAGER_GIGE returns the handle to a new IMAGER_GIGE or the handle to
%      the existing singleton*.
%
%      IMAGER_GIGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGER_GIGE.M with the given input arguments.
%
%      IMAGER_GIGE('Property','Value',...) creates a new IMAGER_GIGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imager_GigE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imager_GigE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imager_GigE

% Last Modified by GUIDE v2.5 01-Aug-2018 09:31:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imager_GigE_OpeningFcn, ...
                   'gui_OutputFcn',  @imager_GigE_OutputFcn, ...
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


% --- Executes just before imager_GigE is made visible.
function imager_GigE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imager_GigE (see VARARGIN)

% Choose default command line output for imager_GigE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imager_GigE wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global scl ACQ LUMserial imager_handle;

imager_handle = handles;

rmpath('C:\ISI acquisition and analysis\imager 2.0')

%%%%%%%%%%%%%%Init grab/save parameters%%%%%%%

ACQ.FPS = 10;  % frames per second
ACQ.bin = 2;
ACQ.timecourseBit = 0; 
ACQ.btwTrialShutter = 0;
ACQ.MechbtwTrialShutter = 0;

ACQ.servoCloseVal = 0.55;
ACQ.servoOpenVal = 0.3;

ACQ.ROIcrop = [];  %When empty, it saves the full image

ACQ.offLineSpatialBinning = 1;
ACQ.offLineTemporalBinning = 1;

ACQ.sensorGain = 0;
ACQ.Gamma = 1;

ACQ.camera = 'AVT';

% ACQ.x_Offset = 0;
% ACQ.y_Offset = 0;
% ACQ.Xwidth = 0;
% ACQ.Ywidth = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.binning,'value',ACQ.bin);
set(handles.framerate,'String',num2str(ACQ.FPS));

set(handles.offLineSpatialBinning,'String',num2str(ACQ.offLineSpatialBinning));
set(handles.offLineTemporalBinning,'String',num2str(ACQ.offLineTemporalBinning));


%% AVT communication (image acq toolbox)

imaqreset %Clean house
delete(imaqfind) %Clean house

imaqmex('feature','-limitPhysicalMemoryUsage',false);

ACQ.vid = videoinput('gige', 1, 'Mono12');

src = getselectedsource(ACQ.vid);

%Start binning at 1x1 to get chip dimensions
src.BinningHorizontal = 1;
src.BinningVertical = 1;

if strcmp(src.DeviceModelName,'Manta G-235B')
    ACQ.IMGSIZE = [1216 1936];
    ACQ.vid.ROIposition = [0 0 ACQ.IMGSIZE(2) ACQ.IMGSIZE(1)];
end

%Bin it at default above
src.BinningHorizontal = ACQ.bin;
src.BinningVertical = ACQ.bin;

%set frame rate (should be done after the spatial stuff to get accurate reading)
src.AcquisitionFrameRateAbs = ACQ.FPS;
src.ExposureTimeAbs = 1000000/src.AcquisitionFrameRateAbs; %Set exposure to max.  N.B. triggered grabbing will obey maximum possible rate, not 'AcquisitionFrameRateAbs'
dum = round(src.AcquisitionFrameRateAbs*100)/100;%Query actual frame rate
set(handles.framerate,'String',num2str(dum));   %set in GUI
ACQ.FPS = src.AcquisitionFrameRateAbs; %set in structure
%ACQ.FPS = 1/src.ExposureTimeAbs*10^6;  


%'FrameStart' needs to be set prior to setting other trigger properties
src.TriggerSelector = 'FrameStart';  
src.TriggerSource = 'Line1';
src.TriggerActivation = 'LevelHigh';

src.Gain = ACQ.sensorGain;
src.Gamma = ACQ.Gamma;  %Should always be '1'... i.e. no Gamma

%ACQ.vid.ROIposition = [ACQ.xOffset ACQ.yOffset ACQ.xsize ACQ.ysize];

ArmTrigger(1)

updateTempReading

%% Serial communication with luminator

try
    LumCOM = 'COM4';
    try
        fclose(instrfind('Port',LumCOM));
    end
    delete(instrfind('Port',LumCOM))
    
    LUMserial = serial(LumCOM,'Terminator','CR',...
        'DataTerminalReady','on','RequestToSend','on','bytesavailablefcnmode','byte');
    fopen(LUMserial);
    
    LUMsend('tt'); %Enable PC control
    LUMsend('mm'); %Open shutter
    LUMsend('jj'); %Enable extended command set
    
    y = LUMsend('dd') %get intensity level
    
    y = str2num(y);
    set(handles.setlight,'value',y/100); %set slider
catch
    disp('No communication with illuminator. May not be on.')
end

try
    ARDCOM = 'COM5';
    try
        fclose(instrfind('Port',ARDCOM));
    end
    delete(instrfind('Port',ARDCOM))
    
     ACQ.ard = arduino(ARDCOM,'uno')
     ACQ.servo = servo(ACQ.ard,'D9')
     writePosition(ACQ.servo,ACQ.servoCloseVal) 
catch
    disp('No communication with arduino or servo shutter')
end

% --- Outputs from this function are returned to the command line.
function varargout = imager_GigE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in binning.
function binning_Callback(hObject, eventdata, handles)
% hObject    handle to binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns binning contents as cell array
%        contents{get(hObject,'Value')} returns selected item from binning

global ACQ

src = getselectedsource(ACQ.vid);

r = questdlg('Reset binning?  You will lose the ROI','Binning','Yes','No','Yes');
if(strcmp(r,'Yes'))
    
    %% First set to max res and size
    
%     src.BinningHorizontal = 1;
%     src.BinningVertical = 1;
%     
%     ACQ.vid.ROIposition = [0 0 ACQ.IMGSIZE(2) ACQ.IMGSIZE(1)];
    
    %% Set binning
    
    idx = get(hObject,'Value');
    
    src = getselectedsource(ACQ.vid);
    
    binOptions = [1 2 4 8];
    ACQ.bin = binOptions(idx);
    
    src.BinningHorizontal = ACQ.bin;
    src.BinningVertical = ACQ.bin;
    
    
    
    %%
    
    %reset frame rate
    src.AcquisitionFrameRateAbs = ACQ.FPS; %Set frame rate to desired
    src.ExposureTimeAbs = 1000000/src.AcquisitionFrameRateAbs; %Set exposure to max.  N.B. triggered grabbing will obey maximum possible rate, not 'AcquisitionFrameRateAbs'
    dum = round(src.AcquisitionFrameRateAbs*100)/100;%Query actual frame rate
    set(handles.framerate,'String',num2str(dum));   %set in GUI
    ACQ.FPS = src.AcquisitionFrameRateAbs; %set in structure
    
else
    set(hObject,'Value',ACQ.bin);
end
%set(1,'Name','imager_GigE :: PLEASE WAIT! ::'); drawnow;

%triggerconfig(ACQ.vid, 'hardware', 'risingEdge', 'opto0-DB9');

%ACQ.vid.TriggerRepeat = Inf; %Needs to be infinity... stupid im aq tb 



% --- Executes during object creation, after setting all properties.
function binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function framerate_Callback(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerate as text
%        str2double(get(hObject,'String')) returns contents of framerate as a double

global ACQ

ACQ.FPS = str2num(get(hObject,'String'));

src = getselectedsource(ACQ.vid);

src.AcquisitionFrameRateAbs = ACQ.FPS; %Set frame rate to desired

src.ExposureTimeAbs = 1000000/src.AcquisitionFrameRateAbs; %Set exposure to max.  N.B. triggered grabbing will obey maximum possible rate, not 'AcquisitionFrameRateAbs'

dum = round(src.AcquisitionFrameRateAbs*100)/100;%Query actual frame rate
set(handles.framerate,'String',num2str(dum));   %set in GUI
ACQ.FPS = src.AcquisitionFrameRateAbs; %set in structure



% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in focus.
function focus_Callback(hObject, eventdata, handles)
% hObject    handle to focus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ

updateTempReading

RGBflag = get(handles.RGBtag,'value');
if RGBflag
    figure(1)
    c = jet;
    c = interp1(1:64,c,linspace(1,64,256)); %It needs 256 levels... To match the monitor?
    c(end,:) = 1; %Show saturation
    close
end

X = get(handles.focus,'String');

if strcmp(X,'Focus')
    
    %turn off trigger to preview. Its slow to respond otherwise
    ArmTrigger(0)
    
    set(handles.focus,'String','Stop');
    set(handles.Status,'String','Previewing');
    set(handles.Status,'ForegroundColor',[1 0 0]);
    set(handles.contImage,'xtick',[],'ytick',[])
    axes(handles.contImage);     %Make current figure
    himage = image(zeros([ACQ.IMGSIZE/ACQ.bin 3]));
    preview(ACQ.vid,himage);
    if RGBflag
        colormap(c)
    end
elseif strcmp(X,'Stop')
        
    stoppreview(ACQ.vid);
    set(handles.Status,'String','Stopped');
    set(handles.Status,'ForegroundColor',[0 0 0]);
    set(handles.focus,'String','Focus');
    
   
    ArmTrigger(1)
    
end


% --- Executes on button press in grab.
function grab_Callback(hObject, eventdata, handles)
% hObject    handle to grab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ

ACQ.vid.FramesPerTrigger = str2num(get(handles.Ngrab,'string'));
%grab.img = getsnapshot(ACQ.vid);

ArmTrigger(0) %disarm trigger, set to 'immediate'.

start(ACQ.vid)

display('logging data')
while islogging(ACQ.vid)
    
end
display('done logging')

grab.img = getdata(ACQ.vid);
grab.img = squeeze(grab.img);

figure(10);
imagesc(grab.img(:,:,1)),axis off, colormap gray; truesize
r = questdlg('Do you want to save it?','Single Grab','Yes','No','Yes');
if(strcmp(r,'Yes'))
    
    grab.comment = inputdlg('Please enter description:','Image Grab',1,{'No description'},'on');
    
    %These are fields of Stimulator
    animal = get(findobj('Tag','animal'),'String');
    unit   = get(findobj('Tag','unitcb'),'String');
    expt   = get(findobj('Tag','exptcb'),'String');
    datadir= get(findobj('Tag','dataRoot'),'String');
    
    dd = [datadir '\' lower(animal) '\grabs\'];
    if(~exist(dd))
        mkdir(dd);
    end
    fname = [dd 'grab_' animal '_' ...
        unit '_' ...
        expt '_' ...
        datestr(now)];
    fname = strrep(fname,' ','_');
    fname = strrep(fname,':','_');
    fname = strrep(fname,'-','_');
    fname = [fname '.mat'];
    fname(2) = ':';
    %save(fname,'grab');
    
    uisave('grab',fname)
     
end
delete(10);

ArmTrigger(1) %Re-arm hardware trigger

% --- Executes on button press in selectROI.
function selectROI_Callback(hObject, eventdata, handles)
% hObject    handle to selectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ;

ArmTrigger(0)

%Reset to full size before cropping
ACQ.vid.ROIposition = [0 0 ACQ.IMGSIZE(2)/ACQ.bin ACQ.IMGSIZE(1)/ACQ.bin];

img = getsnapshot(ACQ.vid);

figure(10);
imagesc(img), colormap gray; truesize

r = questdlg('Crop the image that is saved to the disk?','Select ROI','Yes','No','Yes');
if(strcmp(r,'Yes'))
    [I2 ROIcrop] = imcrop;
    ACQ.ROIcrop = floor(ROIcrop);
end

ACQ.vid.ROIposition = ACQ.ROIcrop;

close(10)

%save('C:\ISI acquisition and analysis\imager_GigE 2.0\lastROI','ROIcrop')

ArmTrigger(1)

% --- Executes on button press in getLastROI.
function getLastROI_Callback(hObject, eventdata, handles)
% hObject    handle to getLastROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ

[file path] = uigetfile({'*.analyzer'},'Load analyzer file','C:\');

id = find(file == '.');
fext = file(id+1:end);

if file  %if 'cancel' was not pressed
    file = [path file];
        
    load(file,'-mat','Analyzer')
    
    if isfield(Analyzer,'ACQ') 
        ACQ.ROIcrop = Analyzer.ACQ.ROIcrop;    
        display(['Cropping to ' num2str(ACQ.ROIcrop(3)) ' x ' num2str(ACQ.ROIcrop(4))])
        
        ACQ.vid.ROIposition = ACQ.ROIcrop;
    else       
        display('No cropping in selected experiment. Setting to full width')
        ACQ.ROIcrop = [];
    end
    
   
end

% --- Executes on button press in fullROI.
function fullROI_Callback(hObject, eventdata, handles)
% hObject    handle to fullROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ

r = questdlg('Remove cropping window?','No','Yes');
if(strcmp(r,'Yes'))
    ACQ.ROIcrop = [];
    %Reset to full size before cropping
    ACQ.vid.ROIposition = [0 0 ACQ.IMGSIZE(2)/ACQ.bin ACQ.IMGSIZE(1)/ACQ.bin];
end

% --- Executes on slider movement.
function setlight_Callback(hObject, eventdata, handles)
% hObject    handle to setlight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

percPower = round(get(handles.setlight,'value')/255*100);

if percPower<5
    percPower = 5;
end

str = ['d' num2str(percPower)];

LUMsend(str);

% --- Executes during object creation, after setting all properties.
function setlight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setlight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in Shutter.
function Shutter_Callback(hObject, eventdata, handles) %look here
% hObject    handle to Shutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Shutter

ShutterBit = get(handles.Shutter,'value');

if ShutterBit
    str = 'zz\r'; %close
else
    str = 'mm\r'; %open
end

LUMsend(str)

function ISIdataRoot_Callback(hObject, eventdata, handles)
% hObject    handle to ISIdataRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ISIdataRoot as text
%        str2double(get(hObject,'String')) returns contents of ISIdataRoot as a double


% --- Executes during object creation, after setting all properties.
function ISIdataRoot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ISIdataRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function analyzerRoot_Callback(hObject, eventdata, handles)
% hObject    handle to analyzerRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analyzerRoot as text
%        str2double(get(hObject,'String')) returns contents of analyzerRoot as a double


% --- Executes during object creation, after setting all properties.
function analyzerRoot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analyzerRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unittxt_Callback(hObject, eventdata, handles)
% hObject    handle to unittxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unittxt as text
%        str2double(get(hObject,'String')) returns contents of unittxt as a double


% --- Executes during object creation, after setting all properties.
function unittxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unittxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function expttxt_Callback(hObject, eventdata, handles)
% hObject    handle to expttxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expttxt as text
%        str2double(get(hObject,'String')) returns contents of expttxt as a double


% --- Executes during object creation, after setting all properties.
function expttxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expttxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function animaltxt_Callback(hObject, eventdata, handles)
% hObject    handle to animaltxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of animaltxt as text
%        str2double(get(hObject,'String')) returns contents of animaltxt as a double


% --- Executes during object creation, after setting all properties.
function animaltxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to animaltxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in RGBtag.
function RGBtag_Callback(hObject, eventdata, handles)
% hObject    handle to RGBtag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RGBtag


% --- Executes on button press in shutterBetweenTrials.
function shutterBetweenTrials_Callback(hObject, eventdata, handles)
% hObject    handle to shutterBetweenTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shutterBetweenTrials

global ACQ

ACQ.btwTrialShutter = get(handles.shutterBetweenTrials,'value');


% --- Executes on button press in timecourse.
function timecourse_Callback(hObject, eventdata, handles)
% hObject    handle to timecourse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timecourse

global ACQ

ACQ.timecourseBit = get(handles.timecourse,'value'); 

if ACQ.timecourseBit
    
    ArmTrigger(0)
    img = getsnapshot(ACQ.vid);
    
    figure(10);
    imagesc(img),axis off, colormap gray; truesize
    
    r = questdlg('Create ROI for intertrial timecourse analysis?','Select ROI','Yes','No','Yes');
    if(strcmp(r,'Yes'))
        ACQ.timecourseROI = roipoly;
    else
        set(handles.timecourse,'value',0);
        ACQ.timecourseBit = 0;
    end
    ArmTrigger(1)
    close(10)
end

%ResetMatrox

function y = clsend(str)

%Send stuff to Dalsa over camera link serial comm
global scl

%scl = instrfind('Tag','clser');
if(~isempty(scl))
    fprintf(scl,'%s\n',str);
    pause(0.05);
    N = get(scl,'BytesAvailable');
    y = [];
    while(N>0)
        y = [y char(fread(scl,N,'char')')];
        pause(0.05);
        N = get(scl,'BytesAvailable');
    end
else
    y = '\n Error>> No message from camera!\n';
end



function y = LUMsend(str)

%Send stuff to the luminator
global LUMserial

if(~isempty(LUMserial))
    fwrite(LUMserial,[str 13]);
    pause(0.2);
    N = get(LUMserial,'BytesAvailable');
    if(N>0)
        y = fgetl(LUMserial);
    end
else
    y = '\n Error>> No message from luminator!\n';
end


function ArmTrigger(flag)

global ACQ

src = getselectedsource(ACQ.vid);
stop(ACQ.vid)
if ~flag
    triggerconfig(ACQ.vid, 'immediate');
    src.TriggerMode = 'Off';    
     ACQ.vid.TriggerRepeat = 0;
else    
    triggerconfig(ACQ.vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
    %triggerconfig(ACQ.vid, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
    src.TriggerMode = 'On';    
    ACQ.vid.TriggerRepeat = inf;
end


% --- Executes on button press in mechShutter.
function mechShutter_Callback(hObject, eventdata, handles)
% hObject    handle to mechShutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mechShutter

global ACQ

ShutterBit = get(handles.mechShutter,'value');

if ShutterBit
    writePosition(ACQ.servo,ACQ.servoCloseVal) %closed
else
    writePosition(ACQ.servo,ACQ.servoOpenVal) %open
end



% --- Executes on button press in mechShutterBetweenTrials.
function mechShutterBetweenTrials_Callback(hObject, eventdata, handles)
% hObject    handle to mechShutterBetweenTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mechShutterBetweenTrials

global ACQ

ACQ.MechbtwTrialShutter = get(handles.mechShutterBetweenTrials,'value');



function Ngrab_Callback(hObject, eventdata, handles)
% hObject    handle to Ngrab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ngrab as text
%        str2double(get(hObject,'String')) returns contents of Ngrab as a double


% --- Executes during object creation, after setting all properties.
function Ngrab_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ngrab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offLineSpatialBinning_Callback(hObject, eventdata, handles)
% hObject    handle to offLineSpatialBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offLineSpatialBinning as text
%        str2double(get(hObject,'String')) returns contents of offLineSpatialBinning as a double

global ACQ

ACQ.offLineSpatialBinning = str2num(get(handles.offLineSpatialBinning,'String'));

% --- Executes during object creation, after setting all properties.
function offLineSpatialBinning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offLineSpatialBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offLineTemporalBinning_Callback(hObject, eventdata, handles)
% hObject    handle to offLineTemporalBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offLineTemporalBinning as text
%        str2double(get(hObject,'String')) returns contents of offLineTemporalBinning as a double

global ACQ

ACQ.offLineTemporalBinning = str2num(get(handles.offLineTemporalBinning,'String'));


% --- Executes during object creation, after setting all properties.
function offLineTemporalBinning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offLineTemporalBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in RGBtag.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to RGBtag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RGBtag


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on slider movement.
function sensorGain_Callback(hObject, eventdata, handles)
% hObject    handle to sensorGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global ACQ

get(handles.sensorGain,'value')

ACQ.sensorGain = round(get(handles.sensorGain,'value')*40);

src = getselectedsource(ACQ.vid);
src.Gain = ACQ.sensorGain;

% --- Executes during object creation, after setting all properties.
function sensorGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensorGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
