function varargout = MW(varargin)
% MW MATLAB code for MW.fig
%      MW, by itself, creates a new MW or raises the existing
%      singleton*.
%
%      H = MW returns the handle to a new MW or the handle to
%      the existing singleton*.
%
%      MW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MW.M with the given input arguments.
%
%      MW('Property','Value',...) creates a new MW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MW_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MW_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MW

% Last Modified by GUIDE v2.5 17-May-2017 10:43:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MW_OpeningFcn, ...
                   'gui_OutputFcn',  @MW_OutputFcn, ...
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


% --- Executes just before MW is made visible.
function MW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MW (see VARARGIN)

% Choose default command line output for MW
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MW wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global GUIhandles Mstate shutterState

Mstate.running = 0;

%Set GUI to the default established in configureMstate
set(handles.screendistance,'string',num2str(Mstate.screenDist))
set(handles.analyzerRoots,'string',Mstate.analyzerRoot)
set(handles.dataRoot,'string',Mstate.dataRoot)
set(handles.animal,'string',Mstate.anim)
set(handles.unitcb,'string',Mstate.unit)
set(handles.exptcb,'string',Mstate.expt)
set(handles.hemisphere,'string',Mstate.hemi)
set(handles.screendistance,'string',Mstate.screenDist)
set(handles.monitor,'string',Mstate.monitor)
set(handles.stimulusIDP,'string',Mstate.stimulusUDP)
set(handles.twophotonflag,'value',Mstate.twoP)
set(handles.widefieldflag,'value',Mstate.WF)
set(handles.ephysflag,'value',Mstate.Ephys)

updateMonitorValues %I don't think this is necessary since its called before building stimuli. It does not send anything over, but just sets the monitor values

GUIhandles.main = handles;

%initialize eye shutter settings
shutterState.use=0;
shutterState.ini=0;


% --- Outputs from this function are returned to the command line.
function varargout = MW_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function animal_Callback(hObject, eventdata, handles)
% hObject    handle to animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of animal as text
%        str2double(get(hObject,'String')) returns contents of animal as a double

global Mstate

Mstate.anim = get(handles.animal,'string');

anaroot = get(handles.analyzerRoots,'string');
Mstate.analyzerRoot = anaroot;

roots = parseString(Mstate.analyzerRoot,';');

dirinfo = dir([roots{1} '\' Mstate.anim]); %Use the first root path for the logic below

if length(dirinfo) > 2 %If the animal folder exists and there are files in it
    
    lastunit = dirinfo(end).name(6:8);
    lastexpt = dirinfo(end).name(10:12);

    newunit = lastunit; 
    newexpt = sprintf('%03d',str2num(lastexpt)+1); %Go to next experiment number
    
else  %if animal folder does not exist or there aren't any files.  The new folder will
        %be created when you hit the 'run' button
    
    newunit = '000';
    newexpt = '000';

end

Mstate.unit = newunit;
Mstate.expt = newexpt;
set(handles.exptcb,'string',newexpt)
set(handles.unitcb,'string',newunit)

UpdateACQExptName   %Send expt info to acquisition


% --- Executes during object creation, after setting all properties.
function animal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to animal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hemisphere_Callback(hObject, eventdata, handles)
% hObject    handle to hemisphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hemisphere as text
%        str2double(get(hObject,'String')) returns contents of hemisphere as a double

global Mstate

%This is not actually necessary since updateMstate is always called prior
%to showing stimuli...
Mstate.hemi = get(handles.hemisphere,'string');

% --- Executes during object creation, after setting all properties.
function hemisphere_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hemisphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function screendistance_Callback(hObject, eventdata, handles)
% hObject    handle to screendistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of screendistance as text
%        str2double(get(hObject,'String')) returns contents of screendistance as a double

global Mstate

%This is not actually necessary since updateMstate is always called prior
%to showing stimuli...  
Mstate.screenDist = str2num(get(handles.screendistance,'string'));


% --- Executes during object creation, after setting all properties.
function screendistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to screendistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Mstate GUIhandles Pstate trialno analogIN ACQ TwoPcomState RipStruct DcomState

DcomState.serialPortHandle.bytesavailablefcn = @Displaycb; 

%Run it!
if ~Mstate.running
    
    if Mstate.WF
        DcomState.serialPortHandle.bytesavailablefcn = ''; %Not sure this is necessary
        if isrunning(ACQ.vid)
            disp('Cancelling previous acquisition...')
            stop(ACQ.vid)
        end
    end
        
    %Check if this analyzer file already exists!
    roots = parseString(Mstate.analyzerRoot,';');    
    for i = 1:length(roots)  %loop through each root
        title = [Mstate.anim '_' sprintf('u%s',Mstate.unit) '_' Mstate.expt];
        dd = [roots{i} '\' Mstate.anim '\' title '.analyzer'];
        
        if(exist(dd))
            warndlg('Directory exists!!!  Please advance experiment before running')
            return
        end
    end
    
    %Check if this data file already exists!
    if Mstate.WF
        UpdateACQExptName %not necessary, but whatever
        
        [Oflag dd] = checkforOverwrite;
        if Oflag
            return
        else
            dd
            mkdir(dd)
        end
        
    end
    
    Mstate.running = 1;  %Global flag for interrupt in real-time loop ('Abort')
    
    %Update states just in case user has not pressed enter after inputing
    %fields:
    updateLstate
    updateMstate    
    
    makeLoop;  %makes 'looperInfo'.  This must be done before saving the analyzer file.
    
%     if strcmp('RD',getmoduleID) %if raindropper
%         if getParamVal('randseed_bit',1) %if random seed bit is set
%             rval = round(rand(1)*1000)/1000;
%             updatePstate('rseed',rval)
%         end
%     end

    saveExptParams  %Save .analyzer. Do this before running... in case something crashes

    set(handles.runbutton,'string','Abort')    
        
    %%%Get the 2photon Acquisition running:
    if Mstate.twoP;  %Flag for the link with scanimage
        fprintf(TwoPcomState.serialPortHandle,'G');  %GO!!!
        pause(1.5) %Need to give it time to ramp up the res mirror
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%Send initial parameters to display
    sendPinfo
    waitforDisplayResp
    sendMinfo
    waitforDisplayResp
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Set camera parameters
    P = getParamStruct(1);
    if Mstate.WF
        ACQ.total_time = P.predelay+P.postdelay+P.stim_time;

        %This is a hack to account for the fact that frames are only
        %grabbed during "logical high".  Sometimes its shorter than
        %expected on long trials.  Rounding error?
%         if ACQ.total_time>10
%             ACQ.vid.FramesPerTrigger = ceil(ACQ.total_time*ACQ.FPS)-4;  %frames per trial
%         else
%             ACQ.vid.FramesPerTrigger = ceil(ACQ.total_time*ACQ.FPS);
%         end
        
        ACQ.vid.FramesPerTrigger = ceil(ACQ.total_time*ACQ.FPS)-4; 
        
        %imaqmem(10^10);  %Setting it bigger than the limit will set it to the max
        
        %Make sure analog in is not running
%         stop(analogIN)
%         flushdata(analogIN)        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
          
    %%%% for ephys exp %%%%
    
    if Mstate.Ephys
        sendFileNametoRipple; %This only creates a folder
        %         opers = xippmex('opers');
        %         desc = xippmex('trial', opers(1));
        disp('Connecting to Trellis...')
        status = xippmex;
        disp('done')
        
        %[trial_descriptor] = xippmex('trial', RipStruct.opers(1), 'recording', RipStruct.Filename,[],RipStruct.auto_incr,[])
        [trial_descriptor] = xippmex('trial','recording', RipStruct.Filename,[],RipStruct.auto_incr,[])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    trialno = 1;
    
    %In 2 computer version 'run2' is no longer a loop, but gets recalled
    %after each trial... 
    %In 'endAcquisition' (2ph), or 'Displaycb' (intrinsic)

    if Mstate.WF %beginning of first trial
        start(ACQ.vid) %Arm the hardware trigger
    end
    
   
    set(0,'RecursionLimit',1300)
    
    if Mstate.WF
        %run2_beta2
        run2
    else
        run2
    end 
    
    
    
    if Mstate.WF %end of last trial
        disp('Disarming trigger')
        stop(ACQ.vid)
    end
    
    
    
else % 'Abort'
    Mstate.running = 0;  %Global flag for interrupt in real-time loop ('Abort')    
    set(handles.runbutton,'string','Run')    
    
    if Mstate.WF
        disp('Disarming trigger')
        stop(ACQ.vid)
    end
    
    if Mstate.twoP  %Flag for the link with Scanbox
        fprintf(TwoPcomState.serialPortHandle,'S');  %Stop acquiring
    end
    
    if Mstate.Ephys %get(handles.ephysTag,'value');  %Flag for the link with scanimage
        %xippmex('trial', RipStruct.opers(1), 'stopped', RipStruct.Filename,[],RipStruct.auto_incr,[]);
        status = xippmex;
        %xippmex('trial', RipStruct.opers(1), 'stopped')
        xippmex('trial','stopped')
        disp('Abort: Ephys');
        %xippmex('close');
    end
end

%This is done to ensure that user builds a new stimulus before doing
%'playsample'.  Otherwise it will open the shutter.
%set(GUIhandles.param.playSample,'enable','off')


% --- Executes on button press in unitcb.
function unitcb_Callback(hObject, eventdata, handles)
% hObject    handle to unitcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Mstate

newunit = sprintf('%03d',str2num(Mstate.unit)+1);
Mstate.unit = newunit;
set(handles.unitcb,'string',newunit)

newexpt = '000';
Mstate.expt = newexpt;
set(handles.exptcb,'string',newexpt)

UpdateACQExptName   %Send expt info to acquisition

% --- Executes on button press in exptcb.
function exptcb_Callback(hObject, eventdata, handles)
% hObject    handle to exptcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Mstate

newexpt = sprintf('%03d',str2num(Mstate.expt)+1);
Mstate.expt = newexpt;
set(handles.exptcb,'string',newexpt)

UpdateACQExptName   %Send expt info to acquisition


% --- Executes on button press in closeDisplay.
function closeDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to closeDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DcomState

fwrite(DcomState.serialPortHandle,'C;~')


% --- Executes on button press in twophotonflag.
function twophotonflag_Callback(hObject, eventdata, handles)
% hObject    handle to twophotonflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twophotonflag

global GUIhandles Mstate

Mstate.twoP = get(handles.twophotonflag,'value');
set(GUIhandles.main.twophotonflag,'value',Mstate.twoP)



function analyzerRoots_Callback(hObject, eventdata, handles)
% hObject    handle to analyzerRoots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analyzerRoots as text
%        str2double(get(hObject,'String')) returns contents of analyzerRoots as a double


%This is not actually necessary since updateMstate is always called prior
%to showing stimuli...

global Mstate

Mstate.analyzerRoot = get(handles.analyzerRoots,'string');
UpdateACQExptName


% --- Executes during object creation, after setting all properties.
function analyzerRoots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analyzerRoots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in REflag.
function REflag_Callback(hObject, eventdata, handles)
% hObject    handle to REflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of REflag

REbit = get(handles.REflag,'value');
moveShutter(2,REbit)
waitforDisplayResp



% --- Executes on button press in LEflag.
function LEflag_Callback(hObject, eventdata, handles)
% hObject    handle to LEflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LEflag

LEbit = get(handles.LEflag,'value');
moveShutter(1,LEbit)
waitforDisplayResp



function monitor_Callback(hObject, eventdata, handles)
% hObject    handle to monitor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of monitor as text
%        str2double(get(hObject,'String')) returns contents of monitor as a double

global Mstate

Mstate.monitor = get(handles.monitor,'string');

updateMonitorValues
sendMonitor

% --- Executes during object creation, after setting all properties.
function monitor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to monitor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stimulusIDP_Callback(hObject, eventdata, handles)
% hObject    handle to stimulusIDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stimulusIDP as text
%        str2double(get(hObject,'String')) returns contents of stimulusIDP as a double


% --- Executes during object creation, after setting all properties.
function stimulusIDP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimulusIDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in widefieldflag.
function widefieldflag_Callback(hObject, eventdata, handles)
% hObject    handle to widefieldflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of widefieldflag

global GUIhandles Mstate

Mstate.WF = get(handles.widefieldflag,'value');
set(GUIhandles.main.widefieldflag,'value',Mstate.WF)


% --- Executes on button press in F1analysisFlag.
function F1analysisFlag_Callback(hObject, eventdata, handles)
% hObject    handle to F1analysisFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F1analysisFlag

global GUIhandles

flag = get(handles.F1analysisFlag,'value');
set(GUIhandles.main.F1analysisFlag,'value',flag)



% --- Executes on button press in F0analysisFlag.
function F0analysisFlag_Callback(hObject, eventdata, handles)
% hObject    handle to F0analysisFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F0analysisFlag



function dataRoot_Callback(hObject, eventdata, handles)
% hObject    handle to dataRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataRoot as text
%        str2double(get(hObject,'String')) returns contents of dataRoot as a double

global Mstate

Mstate.dataRoot = handles.dataRoot.String;
UpdateACQExptName


% --- Executes during object creation, after setting all properties.
function dataRoot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataRoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ephysflag.
function ephysflag_Callback(hObject, eventdata, handles)
% hObject    handle to ephysflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ephysflag

global GUIhandles Mstate RipStruct


Mstate.Ephys = get(handles.ephysflag,'value');
%set(GUIhandles.main.ephysflag,'value',flag)

% ephysFlag = get(handles.ephysTag,'value');
if Mstate.Ephys %ephysFlag
    
    msgbox('Set the electrode in GUI.');
    msgbox({'See command window to make sure Trellis connected.' 'May need to restart Trellis.'})
    
    path(path,'C:\Program Files (x86)\Ripple\Trellis\Tools\xippmex')
%     xippmex('close');
    % Initialize xippmex
    disp('Connecting to Trellis...')
    status = xippmex;
    disp('done')
    if status ~= 1; 
        error('Xippmex Did Not Initialize');  
    end
    
    RipStruct.elecs = xippmex('elec', 'all');
    [RipStruct.state] = xippmex('signal',RipStruct.elecs);
    
%     opers = xippmex('opers');
%     if ~isempty(opers)
%         desc = xippmex('trial', opers(1));
%        
%         RipStruct.opers = opers; % xippmex('opers'); % Trellis operator ID
%         RipStruct.desc = desc; % xippmex('trial', RipStruct.opers);
%     else
%         disp('error in opers')
%     end
end

function sendFileNametoRipple %(handles,Mstate)

global Mstate RipStruct

%% Send file name to Ripple

%ephysBit = get(handles.ephysTag,'value');
msg = [Mstate.anim '_u' Mstate.unit '_' Mstate.expt];

%RipStruct.Filename = sprintf('C:/ephys/dataFiles/%s',msg);
anim_dir = sprintf('%s/%s',Mstate.dataRoot,Mstate.anim);
mkdir(anim_dir)
RipStruct.Filename = sprintf('%s/%s',anim_dir,msg);

RipStruct.auto_incr = 0;


% --- Executes on selection change in electrodeType.
function electrodeType_Callback(hObject, eventdata, handles)
% hObject    handle to electrodeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns electrodeType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from electrodeType

global Mstate

eflag = get(handles.electrodeType,'value');
switch eflag
    
    case 1
        Mstate.electrode = 'A1x16';
    case 2
        Mstate.electrode = 'A2x16';
    case 3
        Mstate.electrode = 'A4x8';
    case 4
        Mstate.electrode = 'A8x8';
end


% --- Executes during object creation, after setting all properties.
function electrodeType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrodeType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cameraType.
function cameraType_Callback(hObject, eventdata, handles)
% hObject    handle to cameraType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cameraType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cameraType


global Mstate

cflag = get(handles.cameraType,'value');
switch cflag
    
    case 1
        Mstate.camera = 'Dalsa';
    case 2
        Mstate.camera = 'AVT';

end

% --- Executes during object creation, after setting all properties.
function cameraType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cameraType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
