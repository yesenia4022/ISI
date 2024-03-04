function varargout = pF0_ISI(varargin)

%Ian Nauhaus

% PF0_ISI M-file for pF0_ISI.fig
%      PF0_ISI, by itself, creates a new PF0_ISI or raises the existing
%      singleton*.
%
%      H = PF0_ISI returns the handle to a new PF0_ISI or the handle to
%      the existing singleton*.
%
%      PF0_ISI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PF0_ISI.M with the given input arguments.
%
%      PF0_ISI('Property','Value',...) creates a new PF0_ISI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before processF0_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pF0_ISI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pF0_ISI

% Last Modified by GUIDE v2.5 28-Jun-2017 14:51:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pF0_ISI_OpeningFcn, ...
                   'gui_OutputFcn',  @pF0_ISI_OutputFcn, ...
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


% --- Executes just before pF0_ISI is made visible.
function pF0_ISI_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pF0_ISI (see VARARGIN)

%Folders for the "other" processF0
% rmpath('F:\neurostuff\2phAnalysis_pep\AnalysisCode')
% rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI')
% rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI\New')
% rmpath('F:\neurostuff\2phAnalysis_pep\cbpep')
% rmpath('F:\neurostuff\2phAnalysis_pep\pep2matlab_new')
% rmpath('F:\neurostuff\2phAnalysis_pep\pepanalysis')
% 
% 
% %Folders for the "other" processF0
% rmpath('F:\neurostuff\2phAnalysis\AnalysisCode')
% rmpath('F:\neurostuff\2phAnalysis\2ph_Processing')
% rmpath('F:\neurostuff\2phAnalysis\AnalysisCode\DynamicProcess')
% rmpath('F:\neurostuff\2phAnalysis\2pAnGUI')
% rmpath('F:\neurostuff\2phAnalysis\2pAnGUI\general')
% 
% %Folders for the "other" pF0_ISI
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\2pAnGUI')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\2pAnGUI\New')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\cbpep')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\pep2matlab_new')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\pepanalysis')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode\ContrastResp')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode\Masking')
% rmpath('C:\2ph_code\Beta\2phAnalysis_pep\AnalysisCode\DynamicProcess')

%Folders for this pF0_ISI
p = 'C:\Users\Ian Nauhaus\Documents\Matlab\ISI acquisition and analysis\Processing\';
path([p 'AnalysisCode'],path)
path([p 'ISI_Processing'],path)
path([p 'ISIAnGUI'],path)
path([p 'ISIAnGUI\general'],path)
path([p 'CAnalysisCode\ContrastResp'],path)
path([p 'AnalysisCode\DynamicProcess'],path)
path([p 'AnalysisCode\DynamicProcess\RevCorr_GUI'],path)
path([p 'offlineMovementCorrection'],path)
path([p 'getMouseAreas'],path)
path([p '\Model Fitting'],path)
path([p 'Color vs Ret'],path)


global G_handles

G_handles = handles;

% Choose default command line output for pF0_ISI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pF0_ISI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


clear global f0m f0m_var Tens Tens_var repdom funcmap bcond bsflag maskS symbolInfo 

% --- Outputs from this function are returned to the command line.
function varargout = pF0_ISI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)

%First make the space-time tensors, Tens and Tens_var
%processButton

global Tens Tens_var bsflag repDom G_handles AUE

t0 = cputime;

bsflag = get(G_handles.basesub,'Value');

b = str2double(get(G_handles.bstart,'String')); %in msec as well
b(2) = str2double(get(G_handles.bstop,'String'));

slowMo = get(G_handles.slowMotionFlag,'Value');  %lateral movement correction

repdum = get(G_handles.repDom,'string');
if strcmp(repdum,'All')
    repDom = 1:getnorepeats(1);
else
    eval(['repDom = [' repdum '];'])
end

setacqinfo(1)  %need to set to first trial so that it doesn't use a blank trial when doing CondTensor

set(G_handles.status,'string','Processing...'), drawnow

[Tens Tens_var] = CondTensor3(b,slowMo);  %%Compute entire space time block for each condition

set(G_handles.status,'string','Done'), drawnow
%sound(.6*sin(2*pi*400/(1200)*(0:400)),1200)  %Signal done

t1 = cputime-t0;

set(G_handles.time,'string',num2str(t1))

set(G_handles.loaded,'string',AUE)

set(G_handles.trialMean,'enable','on')

% --- Executes on button press in trialMean.
function trialMean_Callback(hObject, eventdata, handles)
% hObject    handle to trialMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Compute trial mean from Tensor images

global f0m f0m_var Tens Tens_var Flim G_handles

Flim = str2double(get(G_handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(G_handles.epistop,'String')); %Frame stop in ms (to average)

if ~get(G_handles.F1flag,'Value')
    f0m = CondF0(Tens,Flim); %%%Compute |F1| over time interval [Flim(1) Flim(2)]%%%
    if ~isempty(Tens_var{1})
        f0m_var = CondF0(Tens_var,Flim); %%%Compute |F1| over time interval [Flim(1) Flim(2)]%%%
    end
else
    f0m = CondF1(Tens,Flim); %%%Compute mean over time interval [Flim(1) Flim(2)]%%%
    if ~isempty(Tens_var{1})
        f0m_var = CondF1(Tens_var,Flim); %%%Compute |F1| over time interval [Flim(1) Flim(2)]%%%
    end
end

set(G_handles.plot,'enable','on')

set(G_handles.save,'enable','on')

% --- Executes on selection change in func.
function func_Callback(hObject, eventdata, handles)
% hObject    handle to func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns func contents as cell array
%        contents{get(hObject,'Value')} returns selected item from func


% --- Executes during object creation, after setting all properties.
function func_CreateFcn(hObject, eventdata, handles)
% hObject    handle to func (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
global f0m Tens cellS
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');

if filename ~= 0
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    
    if isfield(S,'f0m')
        f0m = S.f0m;    %f0m is a cell array with images from each condition
        Tens = S.Tens;
    end
    cellS = S.cellS;
    set(handles.plot,'enable','on')
    
    set(handles.loaded,'string',filename(1:length(filename)-4))

end

function setimagedir_Callback(hObject, eventdata, handles)
% hObject    handle to setimagedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setimagedir as text
%        str2double(get(hObject,'String')) returns contents of setimagedir as a double


function loadexp_Callback(hObject, eventdata, handles)
% hObject    handle to loadexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadexp as text
%        str2double(get(hObject,'String')) returns contents of loadexp as a double


% --- Executes during object creation, after setting all properties.
function loadexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setdirs.
function setdirs_Callback(hObject, eventdata, handles)
% hObject    handle to setdirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ACQinfo Analyzer maskS G_handles

Gsetdirectories

fno = 1; tno = 1;
set(handles.frameno,'string',num2str(fno))
set(handles.trialno,'string',num2str(tno))
[Im] = Load2phImage(fno,tno);

axes(G_handles.rimage1);     %Make rimage current figure
cla
imagesc(Im), colormap gray 
set(handles.rimage1,'xtick',[],'ytick',[])

setGUIlabels

%handle the slider
nF = ACQinfo.numberOfFrames;
nT = getnotrials;
set(G_handles.frameslide,'value',(fno-1)/(nF-1))
set(G_handles.trialslide,'value',(tno-1)/(nT-1))

stepsize = 1/nF;
set(G_handles.frameslide,'SliderStep',[stepsize 4*stepsize])
stepsize = 1/nT;
set(G_handles.trialslide,'SliderStep',[stepsize 4*stepsize])


function datadir_Callback(hObject, eventdata, handles)
% hObject    handle to datadir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datadir as text
%        str2double(get(hObject,'String')) returns contents of datadir as a double


% --- Executes during object creation, after setting all properties.
function datadir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datadir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epistart_Callback(hObject, eventdata, handles)
% hObject    handle to epistart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epistart as text
%        str2double(get(hObject,'String')) returns contents of epistart as a double


% --- Executes during object creation, after setting all properties.
function epistart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epistart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epistop_Callback(hObject, eventdata, handles)
% hObject    handle to epistop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epistop as text
%        str2double(get(hObject,'String')) returns contents of epistop as a double


% --- Executes during object creation, after setting all properties.
function epistop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epistop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_Callback(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau as text
%        str2double(get(hObject,'String')) returns contents of tau as a double


% --- Executes during object creation, after setting all properties.
function tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HPBW_Callback(hObject, eventdata, handles)
% hObject    handle to HPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HPBW as text
%        str2double(get(hObject,'String')) returns contents of HPBW as a double


% --- Executes during object creation, after setting all properties.
function HPBW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LPBW_Callback(hObject, eventdata, handles)
% hObject    handle to LPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LPBW as text
%        str2double(get(hObject,'String')) returns contents of LPBW as a double


% --- Executes during object creation, after setting all properties.
function LPBW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LPBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bstart_Callback(hObject, eventdata, handles)
% hObject    handle to bstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bstart as text
%        str2double(get(hObject,'String')) returns contents of bstart as a double


% --- Executes during object creation, after setting all properties.
function bstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bstop_Callback(hObject, eventdata, handles)
% hObject    handle to bstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bstop as text
%        str2double(get(hObject,'String')) returns contents of bstop as a double


% --- Executes during object creation, after setting all properties.
function bstop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in basesub.
function basesub_Callback(hObject, eventdata, handles)
% hObject    handle to basesub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of basesub

% --- Executes on button press in tempfilt.
function tempfilt_Callback(hObject, eventdata, handles)
% hObject    handle to tempfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tempfilt


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


% --- Executes on button press in setROI.
function setROI_Callback(hObject, eventdata, handles)
global bw f0m
% hObject    handle to setROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fno = str2double(get(handles.frameno,'String'));    %Get frame number
tno = str2double(get(handles.trialno,'String'));    %Get frame number

[Im] = Load2phImage(fno,tno);

figure,imagesc(Im), colormap gray        

bw = roipoly;
close

if ~isempty(f0m)
    set(handles.plot,'enable','on')
end


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
global bw f0m funcmap bcond ACQinfo maskS TCWin symbolInfo Analyzer
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(bw) || length(bw(:,1)) ~= ACQinfo.linesPerFrame || length(bw(1,:)) ~= ACQinfo.pixelsPerLine
    bw = ones(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
end

TCWin = str2double(get(handles.TCWin,'string'));

checkvect = [get(handles.F0im,'value') get(handles.funcim,'value')];

hh = makeMapFilter;

anatomyflag = get(handles.anatomyFlag,'value');

setsymbolstruct %Put all the symbol information into global structure

%%Filter raw F0 images with hh and create the functional maps...
if checkvect(2) == 1        %if "Functional Images" is checked
    
%     if cellmaskflag && ~isempty(maskS.bwCell{1}) 
%         [kernPop popBlank CoM pardom] = GetROIKernels(f0m,maskS.bwCell{1});
%     end
    
    switch symbolInfo.domType

        case 1        %domain type is 'axis'
            funcmap = GprocessAxis(f0m,hh);  %output is a vector image

        case 2        %domain type is 'direction'
            funcmap = GprocessDir(f0m,hh);  %output is a vector image

        case 3        %domain type is 'log'
            %funcmap = GprocessLog(f0m,bw,hh);  %output is complex image: real(funcmap) = mag; imag(funcmap) = pref
            funcmap = GprocessLog2(f0m,bw,hh);  %output is complex image: real(funcmap) = mag; imag(funcmap) = pref

        case 4        %domain type is '2D'
            [angx magx angy magy] = Gprocessret(f0m,hh);
            
        case 5        %domain type is 'binary'
            funcmap = GprocessBinary(f0m,hh);
    end

%     if funcflag == 6        %%functionality is color
%         [funcmap colorSel] = GprocessColor(f0m,hh);
%     end
%     
%     if funcflag == 7        %%functionality is color & form
%         [funcmap colorSel] = GprocessColorForm(f0m,hh);
%     end
end


%Create plots

switch symbolInfo.domType

    case 1    %"Axis"

        if checkvect(1) == 1  %mean image for each condition
            plotF0dum(f0m,bw,hh)
        end
        if checkvect(2) == 1  %functional maps
            ang = angle(funcmap);
            mag = abs(funcmap);
            ang = (ang+pi*(1-sign(ang)))/2*180/pi;  %Put domain as [0 180].
           
            Gplotaxismap(bw,real(ang),anatomyflag), title(symbolInfo.str{1},'FontWeight','bold','FontSize',15);

        end

    case 2    %"direction"
        
        if checkvect(1) == 1
            plotF0(f0m,bw,hh)
        end
        if checkvect(2) == 1  %functional maps
            ang = angle(funcmap);
            mag = abs(funcmap);
            ang = (ang+pi*(1-sign(ang)))*180/pi;  %make it 0 to 360.
            
            Gplotdirmap(mag.*bw,ang,anatomyflag), title('Direction','FontWeight','bold','FontSize',15);
                
            %GplotorimapDots(ang,celllocs1)   %To use this other things
            %must change
        end

    case 3    %"log domain"
        
        if checkvect(1) == 1
            plotF0(f0m,bw,hh)
        end
        if checkvect(2) == 1  %functional maps
            
            Int = real(funcmap);
            pref = imag(funcmap);
            
            Gplotlogmap(Int.*bw,pref,anatomyflag), title(symbolInfo.str{1}(find(symbolInfo.str{1}~='_')),'FontWeight','bold','FontSize',15);

        end

    case 4    %Retinotopy
        
        if checkvect(1) == 1
            
            mima = [min(f0m{1}(:)) max(f0m{1}(:))];
            N = length(f0m);
            k = 1;
            figure
            for i = 1:N
                if i ~= bcond+1
                    subplot(1,N-length(bcond),k)
                    imagesc(f0m{i},mima)
                    title(['Condition ' num2str(k-1)])

                    k = k+1;
                end
            end
            colormap gray
        end
        
        if checkvect(2) == 1  %functional maps

            screenDist = Analyzer.M.screenDist;
            screenResX = Analyzer.M.xpixels/Analyzer.M.screenXcm;  %pix/cm
            screenResY = Analyzer.M.ypixels/Analyzer.M.screenYcm;

            [xpos ypos xsize ysize] = getPosSize;
            xsize_cm = (max(xpos)-min(xpos)+min(xsize))/screenResX;  %cm stimulus width
            %xsize_deg = 2*atan2(xsize_cm/2,screenDist)*180/pi;  %convert to deg
            xsize_deg = 360*xsize_cm/(2*pi*screenDist);
            
            ysize_cm = (max(ypos)-min(ypos)+min(ysize))/screenResY;  %cm stimulus width
            %ysize_deg = 2*atan2(ysize_cm/2,screenDist)*180/pi;  %convert to deg
            ysize_deg = 360*ysize_cm/(2*pi*screenDist);

            figure('Name','RETINOTOPIC MAPS','NumberTitle','off')
            subplot(2,1,1), imagesc(angx,'AlphaData',magx.*bw,[0 xsize_deg]); colorbar, set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]), ylabel('x position (deg)')
            subplot(2,1,2), imagesc(angy,'AlphaData',magy.*bw,[0 ysize_deg]); colorbar, set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]), ylabel('y position (deg)')
            plotpixel_cb

            figure('Name','SCREEN COVERAGE','NumberTitle','off')
            Gplotretcoverage(angx,magx,angy,magy)

        end
        
    case 5    %"binary domain" (e.g. ocular dominance)

        if checkvect(1) == 1
            plotF0(f0m,bw,hh)
        end
        
        if checkvect(2) == 1  %functional maps
            
            Gplotbinarymap(funcmap,anatomyflag)
            
            id = find(bw);
            mi = prctile(funcmap(id).*bw(id),2); ma = prctile(funcmap(id).*bw(id),98);
            figure, imagesc(funcmap.*bw,[mi ma]), colormap gray
            title(symbolInfo.str{1}(find(symbolInfo.str{1}~='_')),'FontWeight','bold','FontSize',15);
            axis square
        end

%     case 6    %"color"
%         
%         if checkvect(1) == 1
%             plotF0(f0m,bw,hh)
%         end
%         if checkvect(2) == 1  %functional maps
%             %funcmap = Cellurize(funcmap,maskS.bwCell{1});
%             ang = angle(funcmap);
%             mag = abs(funcmap);
%             ang = (ang+pi*(1-sign(ang)))/2*180/pi;  %Put in orientation domain and convert to degrees.
%             figure
%             Gplotcolormap(bw.*mag,ang), title('Color Axis','FontWeight','bold','FontSize',15);
% 
%             figure
%             subplot(1,2,1)
%             imagesc(colorSel), colormap gray, colorbar
%             colorbar('YTick',[-1 0 1],'YTickLabel',{'-1 Lum','0 Both','+1 Color'})
%             title('Color to Lum Response','FontWeight','bold','FontSize',12);
% 
%         end

end


% --- Executes on button press in F0im.
function F0im_Callback(hObject, eventdata, handles)
% hObject    handle to F0im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F0im



% --- Executes on button press in funcim.
function funcim_Callback(hObject, eventdata, handles)
% hObject    handle to funcim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of funcim



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)

% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global f0m AUE Tens cellS
% 
% UE = get(handles.loadexp,'string');
% path = 'c:\Processed Data\';
% filename = strcat(path,AUE)
% uisave({'f0m','Tens','cellS'},filename)

global f0m AUE Tens

path = 'c:\';
filename = strcat(path,AUE,'_f0m');
uisave('f0m',filename)


function analyzedir_Callback(hObject, eventdata, handles)
% hObject    handle to analyzedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analyzedir as text
%        str2double(get(hObject,'String')) returns contents of analyzedir as a double


% --- Executes during object creation, after setting all properties.
function analyzedir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analyzedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loadana_Callback(hObject, eventdata, handles)
% hObject    handle to loadana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadana as text
%        str2double(get(hObject,'String')) returns contents of loadana as a double


% --- Executes during object creation, after setting all properties.
function loadana_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in lumFlag.
function lumFlag_Callback(hObject, eventdata, handles)
% hObject    handle to lumFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lumFlag


% --- Executes on button press in slowMotionFlag.
function slowMotionFlag_Callback(hObject, eventdata, handles)
% hObject    handle to slowMotionFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slowMotionFlag




% --- Executes on button press in cellmaskflag.
function cellmaskflag_Callback(hObject, eventdata, handles)
% hObject    handle to cellmaskflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cellmaskflag

%at least one has to be checked...
if ~get(handles.cellmaskflag,'value')
    set(handles.F0flag,'value',1)
end



% --- Executes on button press in anatomyFlag.
function anatomyFlag_Callback(hObject, eventdata, handles)
% hObject    handle to anatomyFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anatomyFlag




% --- Executes on button press in cellmaskflag.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to cellmaskflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cellmaskflag


% --- Executes on button press in anatomyFlag.
function radiobutton15_Callback(hObject, eventdata, handles)
% hObject    handle to anatomyFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anatomyFlag



function maskSize_Callback(hObject, eventdata, handles)
% hObject    handle to maskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskSize as text
%        str2double(get(hObject,'String')) returns contents of maskSize as a double


% --- Executes during object creation, after setting all properties.
function maskSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskThresh_Callback(hObject, eventdata, handles)
% hObject    handle to maskThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskThresh as text
%        str2double(get(hObject,'String')) returns contents of maskThresh as a double

resetMask

% --- Executes during object creation, after setting all properties.
function maskThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskMorph_Callback(hObject, eventdata, handles)
% hObject    handle to maskMorph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskMorph as text
%        str2double(get(hObject,'String')) returns contents of maskMorph as a double

resetMask

% --- Executes during object creation, after setting all properties.
function maskMorph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskMorph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in EbarFlag.
function EbarFlag_Callback(hObject, eventdata, handles)
% hObject    handle to EbarFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EbarFlag


function repDom_Callback(hObject, eventdata, handles)
% hObject    handle to repDom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repDom as text
%        str2double(get(hObject,'String')) returns contents of repDom as a double


% --- Executes during object creation, after setting all properties.
function repDom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repDom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showMask.
function showMask_Callback(hObject, eventdata, handles)
% hObject    handle to showMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS f0m

%'f0forMask' uses the entire experiment to compute the template for the mask.
%N.B. Process.raw must have already been run to compute the 'f0m' images...
%and w/o baseline subtraction
%%%%%%
f0forMask = 0;
if f0forMask
    
    Tens = zeros(length(f0m{1}(:,1)),length(f0m{1}(1,:)),length(f0m));
    for k = 1:length(f0m)
        Tens(:,:,k) = f0m{k};
    end
    %%%%%%%

    im = max(Tens,[],3);
    %im = mean(Tens,3);

    [maskS.imZ maskS.im] = Znormalize(str2num(get(handles.maskSize,'string')),im);

else
    [maskS.imZ maskS.im] = Znormalize(str2num(get(handles.maskSize,'string')));
end

if ~isfield(maskS,'bw')
    maskS.bw = ones(size(maskS.imZ{1}));
end

resetMask

button = questdlg('Reset region of interest?','ROI','Yes');
if strcmp(button,'Yes')
    maskS.bw = roipoly;
    resetMask
end

resetMask

set(handles.saveMask,'enable','on')



function searchRange_Callback(hObject, eventdata, handles)
% hObject    handle to searchRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of searchRange as text
%        str2double(get(hObject,'String')) returns contents of searchRange as a double


% --- Executes during object creation, after setting all properties.
function searchRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to searchRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function minCellArea_Callback(hObject, eventdata, handles)
% hObject    handle to minCellArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minCellArea as text
%        str2double(get(hObject,'String')) returns contents of minCellArea as a double

resetMask

% --- Executes during object creation, after setting all properties.
function minCellArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minCellArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function TCWin_Callback(hObject, eventdata, handles)
% hObject    handle to TCWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TCWin as text
%        str2double(get(hObject,'String')) returns contents of TCWin as a double


% --- Executes during object creation, after setting all properties.
function TCWin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TCWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function primSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to primSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of primSymbol as text
%        str2double(get(hObject,'String')) returns contents of primSymbol as a double


% --- Executes during object creation, after setting all properties.
function primSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to primSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in domType.
function domType_Callback(hObject, eventdata, handles)
% hObject    handle to domType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns domType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from domType


% --- Executes during object creation, after setting all properties.
function domType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to domType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in secCollapse.
function secCollapse_Callback(hObject, eventdata, handles)
% hObject    handle to secCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns secCollapse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from secCollapse


% --- Executes during object creation, after setting all properties.
function secCollapse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in secSymbol.
function secSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to secSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns secSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from secSymbol


% --- Executes during object creation, after setting all properties.
function secSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tertCollapse.
function tertCollapse_Callback(hObject, eventdata, handles)
% hObject    handle to tertCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tertCollapse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tertCollapse


% --- Executes during object creation, after setting all properties.
function tertCollapse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tertCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tertSymbol.
function tertSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to tertSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tertSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tertSymbol


% --- Executes during object creation, after setting all properties.
function tertSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tertSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in quatSymbol.
function quatSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to quatSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns quatSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quatSymbol


% --- Executes during object creation, after setting all properties.
function quatSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quatSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in quatCollapse.
function quatCollapse_Callback(hObject, eventdata, handles)
% hObject    handle to quatCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns quatCollapse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quatCollapse


% --- Executes during object creation, after setting all properties.
function quatCollapse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quatCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in quintSymbol.
function quintSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to quintSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns quintSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quintSymbol


% --- Executes during object creation, after setting all properties.
function quintSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quintSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in quintCollapse.
function quintCollapse_Callback(hObject, eventdata, handles)
% hObject    handle to quintCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns quintCollapse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from quintCollapse


% --- Executes during object creation, after setting all properties.
function quintCollapse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quintCollapse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function trialslide_Callback(hObject, eventdata, handles)
% hObject    handle to trialslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Nt = getnotrials;
tno = get(handles.trialslide,'value');
tno = round((Nt-1)*tno)+1;

set(handles.trialno,'string',num2str(tno));    %Get trial number
fno = str2double(get(handles.frameno,'string'));    %Get frame number
[Im] = Load2phImage(fno,tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imagesc(Im), colormap gray       %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])


% --- Executes during object creation, after setting all properties.
function trialslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function frameslide_Callback(hObject, eventdata, handles)
% hObject    handle to frameslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global ACQinfo

Nf = ACQinfo.numberOfFrames;
fno = get(handles.frameslide,'value');
fno = round((Nf-1)*fno)+1;

set(handles.frameno,'string',num2str(fno));    %Get trial number
tno = str2double(get(handles.trialno,'string'));    %Get frame number
[Im] = Load2phImage(fno,tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imagesc(Im), colormap gray        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])


% --- Executes during object creation, after setting all properties.
function frameslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in fastMotionFlag.
function fastMotionFlag_Callback(hObject, eventdata, handles)
% hObject    handle to fastMotionFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fastMotionFlag


% --- Executes on selection change in motionModel.
function motionModel_Callback(hObject, eventdata, handles)
% hObject    handle to motionModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns motionModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from motionModel


% --- Executes during object creation, after setting all properties.
function motionModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to motionModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in revCorrGUI.
function revCorrGUI_Callback(hObject, eventdata, handles)
% hObject    handle to revCorrGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pRev


% --- Executes on button press in F0flag.
function F0flag_Callback(hObject, eventdata, handles)
% hObject    handle to F0flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F0flag

%at least one has to be checked...
if ~get(handles.F0flag,'value')
    set(handles.cellmaskflag,'value',1)
end



% --- Executes on button press in loadSyncs.
function loadSyncs_Callback(hObject, eventdata, handles)
% hObject    handle to loadSyncs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadSyncs





function SelectWin_Callback(hObject, eventdata, handles)
% hObject    handle to SelectWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SelectWin as text
%        str2double(get(hObject,'String')) returns contents of SelectWin as a double


% --- Executes during object creation, after setting all properties.
function SelectWin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function blockSize_Callback(hObject, eventdata, handles)
% hObject    handle to blockSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blockSize as text
%        str2double(get(hObject,'String')) returns contents of blockSize as a double


% --- Executes during object creation, after setting all properties.
function blockSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blockSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);

global G_handles G_RChandles ACQinfo Analyzer maskS f0m cellS bw TCWin symbolInfo bcond funcmap

clear global G_handles G_RChandles ACQinfo Analyzer maskS f0m cellS bw TCWin symbolInfo bcond funcmap



% --- Executes on selection change in indicator.
function indicator_Callback(hObject, eventdata, handles)
% hObject    handle to indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns indicator contents as cell array
%        contents{get(hObject,'Value')} returns selected item from indicator


% --- Executes during object creation, after setting all properties.
function indicator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to indicator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in negativeSignalFlag.
function negativeSignalFlag_Callback(hObject, eventdata, handles)
% hObject    handle to negativeSignalFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of negativeSignalFlag


% --- Executes on button press in F1flag.
function F1flag_Callback(hObject, eventdata, handles)
% hObject    handle to F1flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of F1flag




