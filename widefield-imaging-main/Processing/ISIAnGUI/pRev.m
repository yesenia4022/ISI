function varargout = pRev(varargin)

%Ian Nauhaus

% PREV M-file for pRev.fig
%      PREV, by itself, creates a new PREV or raises the existing
%      singleton*.
%
%      H = PREV returns the handle to a new PREV or the handle to
%      the existing singleton*.
%
%      PREV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREV.M with the given input arguments.
%
%      PREV('Property','Value',...) creates a new PREV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pRev_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pRev_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pRev

% Last Modified by GUIDE v2.5 03-Jul-2011 13:13:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pRev_OpeningFcn, ...
                   'gui_OutputFcn',  @pRev_OutputFcn, ...
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


% --- Executes just before pRev is made visible.
function pRev_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pRev (see VARARGIN)

global G_RChandles

G_RChandles = handles;


% Choose default command line output for pRev
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pRev wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pRev_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Analyzer kernels kernblank kernC

tauN = str2num(get(handles.kernelLength,'string'));

switch Analyzer.P.type
    
    case 'FG'

        Gkernelplots4;
        RF = GMakeRF(kernC,tauN);
        
    case 'RD'

        %Grandposplots3;
        Grandposplots2(kernels,tauN)
        
end

function kernelLength_Callback(hObject, eventdata, handles)
% hObject    handle to kernelLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernelLength as text
%        str2double(get(hObject,'String')) returns contents of kernelLength as a double


% --- Executes during object creation, after setting all properties.
function kernelLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernelLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makeKernels.
function makeKernels_Callback(hObject, eventdata, handles)
% hObject    handle to makeKernels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global Analyzer cellS kernels kernblank G_RChandles

hh = makeTemporalfilter;

trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];

switch Analyzer.P.type

    case 'FG'        
        Ggetrevcorrkernel2(cellS.cellMat,trialdom,hh);
        set(handles.plot,'enable','on')
        set(handles.plotMaps,'enable','on')
        
    case 'RD' 
        [kernels kernblank kernIm] = Ggetrandposkernel2(cellS.cellMat,trialdom,hh);
        set(handles.plot,'enable','on')
        
end



% --- Executes on button press in blankNorm.
function blankNorm_Callback(hObject, eventdata, handles)
% hObject    handle to blankNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blankNorm



function Hwidth_Callback(hObject, eventdata, handles)
% hObject    handle to Hwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Hwidth as text
%        str2double(get(hObject,'String')) returns contents of Hwidth as a double


% --- Executes during object creation, after setting all properties.
function Hwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Hwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lwidth_Callback(hObject, eventdata, handles)
% hObject    handle to Lwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lwidth as text
%        str2double(get(hObject,'String')) returns contents of Lwidth as a double


% --- Executes during object creation, after setting all properties.
function Lwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HPflag.
function HPflag_Callback(hObject, eventdata, handles)
% hObject    handle to HPflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HPflag


% --- Executes on button press in LPflag.
function LPflag_Callback(hObject, eventdata, handles)
% hObject    handle to LPflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LPflag


% --- Executes on selection change in HPWind.
function HPWind_Callback(hObject, eventdata, handles)
% hObject    handle to HPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns HPWind contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HPWind


% --- Executes during object creation, after setting all properties.
function HPWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LPWind.
function LPWind_Callback(hObject, eventdata, handles)
% hObject    handle to LPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LPWind contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LPWind


% --- Executes during object creation, after setting all properties.
function LPWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LPWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getMaps.
function getMaps_Callback(hObject, eventdata, handles)
% hObject    handle to getMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Analyzer cellS kernelsIm G_RChandles

hh = makeTemporalfilter;

trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];

tauN = str2num(get(handles.kernelLength,'string'));

switch Analyzer.P.type

    case 'FG'
        [kernelsIm] = Ggetrevcorrkernel_image(trialdom,hh);
        set(handles.plot,'enable','on')
        
    case 'RD' 
        [kernelsIm kernblankIm] = Ggetrandposkernel_image(trialdom,hh);
        set(handles.plot,'enable','on')
        
end

% --- Executes on button press in kernFromMaps.
function kernFromMaps_Callback(hObject, eventdata, handles)
% hObject    handle to kernFromMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Analyzer

switch Analyzer.P.type

    case 'FG'
        plotRevCorrMap
        plotRevCorrMap_log

        Imtokern

    case 'RD' 
        %plotRandPosMap
        plotRandPosMap_xy2
        
end





% --- Executes on button press in plotMaps.
function plotMaps_Callback(hObject, eventdata, handles)
% hObject    handle to plotMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Gplotaxismap_RC

Gplotlogmap_RC



function dropTrials_Callback(hObject, eventdata, handles)
% hObject    handle to dropTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dropTrials as text
%        str2double(get(hObject,'String')) returns contents of dropTrials as a double


% --- Executes during object creation, after setting all properties.
function dropTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dropTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function logfilePath_Callback(hObject, eventdata, handles)
% hObject    handle to logfilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of logfilePath as text
%        str2double(get(hObject,'String')) returns contents of logfilePath as a double


% --- Executes during object creation, after setting all properties.
function logfilePath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to logfilePath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


