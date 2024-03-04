function varargout = imager(varargin)
% IMAGER MATLAB code for imager.fig
%      IMAGER, by itself, creates a new IMAGER or raises the existing
%      singleton*.
%
%      H = IMAGER returns the handle to a new IMAGER or the handle to
%      the existing singleton*.
%
%      IMAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGER.M with the given input arguments.
%
%      IMAGER('Property','Value',...) creates a new IMAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imager

% Last Modified by GUIDE v2.5 16-Mar-2016 21:04:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imager_OpeningFcn, ...
                   'gui_OutputFcn',  @imager_OutputFcn, ...
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


% --- Executes just before imager is made visible.
function imager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imager (see VARARGIN)

% Choose default command line output for imager
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imager wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global scl ACQ LUMserial;

rmpath('C:\ISI acquisition and analysis\imager GigE')

%%%%%%%%%%%%%%Init grab parameters%%%%%%%
ACQ.FPS = 10;  % frames per second
ACQ.IMGSIZE = 1024; %Dalsa
ACQ.bin = 2;
ACQ.timecourseBit = 0; 
ACQ.btwTrialShutter = 0;
ACQ.MechbtwTrialShutter = 0;

ACQ.servoCloseVal = 0.55;
ACQ.servoOpenVal = 0.3;

ACQ.ROIcrop = [];  %When empty, it saves the full image

ACQ.camera = 'Dalsa';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(handles.binning,'value',ACQ.bin);
set(handles.framerate,'String',num2str(ACQ.FPS));

%%  serial communication with camera...
CamCOM = 'COM3';
try
    fclose(instrfind('Port',CamCOM));
end
delete(instrfind('Port',CamCOM))

%delete(instrfind ('Tag','clser'));  %% just in case it crashed
scl = serial(CamCOM,'Tag','clser','Terminator','CR',...
    'DataTerminalReady','on','RequestToSend','on','bytesavailablefcnmode','byte');
fopen(scl);

y = clsend('sem 7'); %Continuous acquisition mode
y = clsend(sprintf('ssf %d',round(ACQ.FPS)));
y = clsend(['sbm ' num2str(ACQ.bin) ' ' num2str(ACQ.bin)]); %X/Y Binning

% y = clsend('ssg 1 4095'); %Continuous acquisition mode
% y = clsend('ssg 2 4095'); %Continuous acquisition mode
% 
% y = clsend('gm 2'); %Continuous acquisition mode

%% Matrox communication (image acq toolbox)

imaqreset %Clean house
delete(imaqfind) %Clean house

floc = 'C:\ISI acquisition and analysis\imager 2.0\New Matrox DCF\';
fname = ['01M60P_' num2str(ACQ.IMGSIZE/ACQ.bin) 'x' num2str(ACQ.IMGSIZE/ACQ.bin) '_12bit2taps_Bin' num2str(ACQ.bin) 'x' num2str(ACQ.bin) '_Async.dcf'];
ACQ.vid = videoinput('matrox', 1, [floc fname]);

triggerconfig(ACQ.vid, 'hardware', 'risingEdge', 'opto0-DB9'); %Set trigger parameters
ACQ.vid.TriggerRepeat = Inf;

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
%      ACQ.ard = arduino('COM5','uno')
%      ACQ.servo = servo(ACQ.ard,9)
%      writePosition(ACQ.servo,.44) %put it in the middle ('closed')
% above is native code, below is code copied from imager_GigE
%
%
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
function varargout = imager_OutputFcn(hObject, eventdata, handles) 
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

global ACQ;

idx = get(hObject,'Value');

binOptions = [1 2 4 8];
ACQ.bin = binOptions(idx);
M = ACQ.IMGSIZE/ACQ.bin;

%set(1,'Name','imager :: PLEASE WAIT! ::'); drawnow;

floc = 'C:\ISI acquisition and analysis\imager 2.0\New Matrox DCF\';
fname = ['01M60P_' num2str(M) 'x' num2str(M) '_12bit2taps_Bin' num2str(ACQ.bin) 'x' num2str(ACQ.bin) '_Async.dcf'];
ACQ.vid = videoinput('matrox', 1, [floc fname]);

y = clsend(['sbm ' num2str(ACQ.bin) ' ' num2str(ACQ.bin)]); %Send to Dalsa

triggerconfig(ACQ.vid, 'hardware', 'risingEdge', 'opto0-DB9');

ACQ.vid.TriggerRepeat = Inf; %Needs to be infinity... stupid im aq tb 

%Allocate?

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
y = clsend(sprintf('ssf %d',round(ACQ.FPS)));


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

RGBflag = get(handles.RGBtag,'value');
if RGBflag
    figure(99)
    c = jet;
    c = interp1(1:64,c,linspace(1,64,256)); %It needs 256 levels... To match the monitor?
    c(end,:) = 1; %Show saturation
    close
end

X = get(handles.focus,'String');

if strcmp(X,'Focus')
    set(handles.focus,'String','Stop');
    set(handles.Status,'String','Previewing');
    set(handles.Status,'ForegroundColor',[1 0 0]);
    set(handles.contImage,'xtick',[],'ytick',[])
    axes(handles.contImage);     %Make current figure
    himage = image(zeros([ACQ.vid.VideoResolution 3]));
    preview(ACQ.vid,himage);
    if RGBflag
        colormap(c)
    end
elseif strcmp(X,'Stop')
    stoppreview(ACQ.vid);
    set(handles.Status,'String','Stopped');
    set(handles.Status,'ForegroundColor',[0 0 0]);
    set(handles.focus,'String','Focus');
    
    
    %Need to reset the hardware after previewing, otherwise it won't be able to arm the hardware
    %trigger:
    ResetMatrox
    
end

% --- Executes on button press in grab.
function grab_Callback(hObject, eventdata, handles)
% hObject    handle to grab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ


grab.img = getsnapshot(ACQ.vid);

figure;
% imagesc(grab.img,[0 4096]),axis off, colormap gray; truesize
imagesc(grab.img),axis off, colormap gray; truesize
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
%delete(10);


%Need to reset the hardware after previewing, otherwise it won't be able to arm the hardware
%trigger:
ResetMatrox
        

% --- Executes on button press in selectROI.
function selectROI_Callback(hObject, eventdata, handles)
% hObject    handle to selectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQ;

img = getsnapshot(ACQ.vid);

figure(10);
imagesc(img),axis off, colormap gray; truesize

r = questdlg('Crop the image that is saved to the disk?','Select ROI','Yes','No','Yes');
if(strcmp(r,'Yes'))
    [I2 ROIcrop] = imcrop;
    ACQ.ROIcrop = round(ROIcrop);
end
close(10)

%save('C:\ISI acquisition and analysis\imager 2.0\lastROI','ROIcrop')

ResetMatrox

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
function Shutter_Callback(hObject, eventdata, handles)
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
    close(10)
end

ResetMatrox

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


function ResetMatrox

global ACQ

M = ACQ.IMGSIZE/ACQ.bin;

imaqreset
floc = 'C:\ISI acquisition and analysis\imager 2.0\New Matrox DCF\';
fname = ['01M60P_' num2str(M) 'x' num2str(M) '_12bit2taps_Bin' num2str(ACQ.bin) 'x' num2str(ACQ.bin) '_Async.dcf'];
ACQ.vid = videoinput('matrox', 1, [floc fname]);

triggerconfig(ACQ.vid, 'hardware', 'risingEdge', 'opto0-DB9'); %Set trigger parameters
ACQ.vid.TriggerRepeat = Inf;


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
