function varargout = pRevCorr(varargin)
% PREVCORR M-file for pRevCorr.fig
%      PREVCORR, by itself, creates a new PREVCORR or raises the existing
%      singleton*.
%
%      H = PREVCORR returns the handle to a new PREVCORR or the handle to
%      the existing singleton*.
%
%      PREVCORR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREVCORR.M with the given input arguments.
%
%      PREVCORR('Property','Value',...) creates a new PREVCORR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before processF0_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pRevCorr_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pRevCorr

% Last Modified by GUIDE v2.5 26-Sep-2010 17:04:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pRevCorr_OpeningFcn, ...
                   'gui_OutputFcn',  @pRevCorr_OutputFcn, ...
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


% --- Executes just before pRevCorr is made visible.
function pRevCorr_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pRevCorr (see VARARGIN)

%Folders for the "other" processF0
rmpath('F:\neurostuff\2phAnalysis_pep\AnalysisCode')
rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI')
rmpath('F:\neurostuff\2phAnalysis_pep\2pAnGUI\New')
rmpath('F:\neurostuff\2phAnalysis_pep\cbpep')
rmpath('F:\neurostuff\2phAnalysis_pep\pep2matlab_new')
rmpath('F:\neurostuff\2phAnalysis_pep\pepanalysis')


%Folders for the "other" processF0
rmpath('F:\neurostuff\2phAnalysis\AnalysisCode')
rmpath('F:\neurostuff\2phAnalysis\2ph_Processing')
rmpath('F:\neurostuff\2phAnalysis\AnalysisCode\DynamicProcess')
rmpath('F:\neurostuff\2phAnalysis\2pAnGUI')
rmpath('F:\neurostuff\2phAnalysis\2pAnGUI\general')

%Folders for this pRevCorr
path('C:\Beta\2phAnalysis\AnalysisCode',path)
path('C:\Beta\2phAnalysis\2ph_Processing',path)
path('C:\Beta\2phAnalysis\2pAnGUI',path)
path('C:\Beta\2phAnalysis\2pAnGUI\general',path)
path('C:\Beta\2phAnalysis\AnalysisCode\ContrastResp',path)
path('C:\Beta\2phAnalysis\AnalysisCode\DynamicProcess',path)
path('C:\Beta\2phAnalysis\AnalysisCode\Masking',path)

global G_handles

G_handles = handles;

% Choose default command line output for pRevCorr
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pRevCorr wait for user response (see UIRESUME)
% uiwait(handles.figure1);


clear f0m f0m_var Tens Tens_var repdom funcmap bcond bsflag maskS symbolInfo 

% --- Outputs from this function are returned to the command line.
function varargout = pRevCorr_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)
global maskS cellMat AUE ACQinfo synctimes bw

if isempty(bw)
    bw = ones(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
end

NT = getnotrials;

t0 = cputime;

CH = GetTrialData([1 1 0 0],1);

if get(handles.fastMotionFlag,'value')    
    [Px_fast Py_fast] = getTrialMotion3(CH{2});
    temp1 = makeGeoTrx(CH{1},Px_fast,Py_fast);
    temp2 = makeGeoTrx(CH{2},Px_fast,Py_fast);
else
    temp1 = CH{1};
    temp2 = CH{2};
end

temp1 = mean(temp1(:,:,2:end-2),3);  %Template for trial-to-trial motion correction
temp2 = mean(temp2(:,:,2:end-2),3);  %Template for trial-to-trial motion correction

Fs = 1000*ACQinfo.pixelsPerLine/ACQinfo.msPerLine;

for T = 1:NT
    T
    CHs = GetTrialData([1 1 1 0],T);  %load second channel because it usually makes for better motion correction
    
    if get(handles.fastMotionFlag,'value')
        [Px_fast Py_fast] = getTrialMotion3(CHs{2});
        CHs{1} = makeGeoTrx(CHs{1},Px_fast,Py_fast);
    end
    
    if get(handles.slowMotionFlag,'value')
        imdum = mean(CHs{2}(:,:,2:end-2),3);
        [mbest nbest] = getShiftVals(imdum,temp2);  %get the transformation for this trial
        for z = 1:length(CHs{1}(1,1,:))
            CHs{1}(:,:,z) = circshift(CHs{1}(:,:,z),[-mbest -nbest]); %transform
        end
    end
    
    %Get time course for each cell       
    
    %[cellMat{T}] = getGeoTrxTimeCourse(CHs{1},Px_fast,Py_fast,maskS.bwCell1);

    maskS.bwCell1 = maskS.bwCell1.*bw;
    
    cellMat{T} = getCelltimecourse(CHs{1},maskS.bwCell1);
    
    CHs{1} = []; %Helps keep it from running out of memory;
    CHs{2} = [];
    
    %Get Synctimes
    dim = size(CHs{3});
    for i = 1:length(CHs{3}(1,1,:))
        dum = CHs{3}(:,:,i)';
        CHs{3}(:,:,i) = reshape(dum(:),dim(1),dim(2));  %Don't want to create new variable, because of memory issues
    end
    
    synctimes{T} = getFlashGraterSynctimes_CRT(CHs{3}(:),Fs);
    
    id = find(synctimes{T} < getparam('predelay')/2);
    synctimes{T}(id) = [];
    id = find(synctimes{T} > getparam('predelay')+getparam('stim_time')+getparam('postdelay')/2);
    synctimes{T}(id) = [];

    %
end


set(handles.status,'string','Done'), drawnow
%sound(.6*sin(2*pi*400/(1200)*(0:400)),1200)  %Signal done

t1 = cputime-t0;

set(handles.time,'string',num2str(t1))

set(handles.loaded,'string',AUE)

set(handles.plot,'enable','on')
set(handles.makeKernels,'enable','on')
set(handles.save,'enable','on')


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
global f0m Tens
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');

if filename ~= 0
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    
   f0m = S.f0m;    %f0m is a cell array with images from each condition
    Tens = S.Tens;
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
global ACQinfo twophDATADIR AUE bw Analyzer

anim = get(handles.loadana,'String');
expt = get(handles.loadexp,'String');

dir = get(handles.analyzedir,'String'); %partial path for analyzer file
setAnalyzerDirectory([dir anim '\'])

loadAnalyzer(expt)

fno = 1; tno = 1;
set(handles.frameno,'string',num2str(fno))
set(handles.trialno,'string',num2str(tno))

dir = get(handles.datadir,'String'); %partial path for .tiff files 

twophDATADIR = [dir anim '\' expt '\'];  %Make path for .tiff files
AUE = [anim '_' expt]; %loaded unit experiment 'u000_000'

setacqinfo_mat(tno)  %Set the global ACQinfo structure (contains trial info)

[Im] = Load2phImage_mat(fno,[1 1 0 0],tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(Im{1},2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(Im{2},1)       %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])


if ~isempty(bw)
    if length(bw(:,1)) ~= length(Im{1}(:,1)) || length(bw(1,:)) ~= length(Im{1}(1,:))
        bw = ones(size(Im{1}));
    end
else
    bw = ones(size(Im{1}));
end

%handle the slider
nF = ACQinfo.numberOfFrames;
nT = getnotrials;
set(handles.frameslide,'value',(fno-1)/(nF-1))
set(handles.trialslide,'value',(tno-1)/(nT-1))

stepsize = 1/nF;
set(handles.frameslide,'SliderStep',[stepsize 4*stepsize])
stepsize = 1/nT;
set(handles.trialslide,'SliderStep',[stepsize 4*stepsize])

set(handles.showMask,'enable','on')
set(handles.process,'enable','on')

set(handles.setROI,'enable','on')

%%%Patch to account for the wrong trial number provided in the Analyzer for
%%%the blank conditions.  N.B. the stimulus was generated ok, but the
%%%Analyzer file is wrong.
%AnalyzerPatch

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

[Im] = Load2phImage_mat(fno,[1 1 0 0],tno);

figure,imagesc(Im{1}), colormap gray        

bw = roipoly;
close

if ~isempty(f0m)
    set(handles.plot,'enable','on')
end

% --- Executes on button press in makeKernels.
function makeKernels_Callback(hObject, eventdata, handles)
% hObject    handle to makeKernels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ACQinfo Analyzer maskS cellMat synctimes kernels kernblank

togstateHP = get(handles.HPflag,'Value');
togstateLP = get(handles.LPflag,'Value');

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

N = length(cellMat{1}(1,:));

if togstateHP == 1
    Hwidth = str2double(get(handles.Hwidth,'string'));
    Hwidth = round(Hwidth/acqPeriod); %convert to samples
    ind = get(handles.HPWind,'value');

    switch ind
        case 1 %Gaussian
            dom = (1:N)-round(N/2);
            H = exp(-dom.^2/(2*Hwidth^2));
            H = -H/sum(H);
            H(round(N/2)) = 1+H(round(N/2));
        case 2 %Hann
            H = zeros(1,N);
            Hd = hann(Hwidth);
            Hd = -Hd./sum(Hd(:));
            Hd(round(Hwidth/2)) = 1+Hd(round(Hwidth/2));
            H(1:Hwidth) = Hd;
        case 3 %Disc
            H = zeros(1,N);
            Hd = -ones(1,Hwidth)/Hwidth;
            Hd(round(Hwidth/2)) = 1+Hd(round(Hwidth/2));
            H(1:Hwidth) = Hd;
    end
    if togstateLP == 0
        hh = ifft(abs(fft(H)));   %Eliminate phase information
    end
end

if togstateLP == 1
    Lwidth = str2double(get(handles.Lwidth,'string'));
    Lwidth = round(Lwidth/acqPeriod); %convert to samples
    ind = get(handles.LPWind,'value');

    switch ind
        case 1
            dom = (1:N)-round(N/2);
            L = exp(-dom.^2/(2*Lwidth^2));
            L = L/sum(L);
        case 2            
            L = zeros(1,N);
            Ld = hann(Lwidth);
            Ld = Ld./sum(Ld(:));
            L(1:Lwidth) = Ld;            
        case 3            
            L = zeros(1,N);
            Ld = ones(1,Lwidth)/Lwidth;
            L(1:Lwidth) = Ld;            
    end
    if togstateHP == 0
        hh = abs(fft(L(:)))';   %Eliminate phase information
    else
        hh = abs(fft(L(:)).*fft(H(:)))';   %Take mag because phase gives a slight shift.
    end
end

if ~or(togstateLP,togstateHP)
    hh = [];
end


tauN = str2num(get(handles.kernelLength,'string'));

switch Analyzer.P.type

    case 'FG'
        [kernels kernblank] = Ggetrevcorrkernel(cellMat,synctimes,maskS.bwCell1,tauN,hh);
        
    case 'RD' 
        [kernels kernblank] = Ggetrandposkernel(cellMat,synctimes,maskS.bwCell1,tauN,hh);
        
end

% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)

% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if isempty(bw) || length(bw(:,1)) ~= ACQinfo.linesPerFrame || length(bw(1,:)) ~= ACQinfo.pixelsPerLine
%     bw = ones(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
% end

global kernels maskS Analyzer

tauN = str2num(get(handles.kernelLength,'string'));

switch Analyzer.P.type
    
    case 'FG'

        Gkernelplots(kernels,maskS.bwCell1,tauN);
        %RF = GMakeRF(kernels,tauN);
        
    case 'RD'

        Grandposplots(kernels,tauN);
        
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
global f0m AUE Tens
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


UE = get(handles.loadexp,'string');
path = 'c:\Processed Data\';
filename = strcat(path,AUE)
uisave({'f0m','Tens'},filename)


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


% --- Executes on button press in trialVarianceFlag.
function trialVarianceFlag_Callback(hObject, eventdata, handles)
% hObject    handle to trialVarianceFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trialVarianceFlag



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

global maskS

[maskS.imZ maskS.im] = Znormalize(str2num(get(handles.maskSize,'string')));

if ~isfield(maskS,'bw')
    maskS.bw = ones(size(maskS.imZ{1}));
end
    
resetMask

button = questdlg('Reset region of interest?','ROI','Yes');
if strcmp(button,'Yes')
    maskS.bw = roipoly;
    resetMask
end

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
[Im] = Load2phImage_mat(fno,[1 1 0 0],tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(Im{1},2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(Im{2},1)        %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])


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
[Im] = Load2phImage_mat(fno,[1 1 0 0],tno);

axes(handles.rimage1);     %Make rimage current figure
cla
imageRGB(Im{1},2)        %Load and plot frame
set(handles.rimage1,'xtick',[],'ytick',[])

axes(handles.rimage2);     %Make rimage current figure
cla
imageRGB(Im{2},1)        %Load and plot frame
set(handles.rimage2,'xtick',[],'ytick',[])


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





% --- Executes on button press in blankNorm.
function blankNorm_Callback(hObject, eventdata, handles)
% hObject    handle to blankNorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blankNorm




% --- Executes on button press in loadMask.
function loadMask_Callback(hObject, eventdata, handles)
% hObject    handle to loadMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');

if filename ~= 0
    
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    
    maskS = S.maskS;    %f0m is a cell array with images from each condition

    set(handles.maskSize,'string',num2str(maskS.maskSize));
    set(handles.maskThresh,'string',num2str(maskS.maskThresh));
    set(handles.maskMorph,'string',num2str(maskS.maskMorph));
    set(handles.minCellArea,'string',num2str(maskS.minCellArea));
    
    figure(40),
    imagesc(maskS.im{1}), colormap gray
    hold on
    contour(maskS.bwCell1,.5,'r')

end


% --- Executes on button press in saveMask.
function saveMask_Callback(hObject, eventdata, handles)
% hObject    handle to saveMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global maskS

maskS.maskSize = str2num(get(handles.maskSize,'string'));
maskS.maskThresh = str2num(get(handles.maskThresh,'string'));
maskS.maskMorph = str2num(get(handles.maskMorph,'string'));
maskS.minCellArea = str2num(get(handles.minCellArea,'string'));

AUE = get(handles.loadexp,'string');
path = 'c:\mask_data\';
filename = strcat(path,AUE)
uisave({'maskS'},filename)


