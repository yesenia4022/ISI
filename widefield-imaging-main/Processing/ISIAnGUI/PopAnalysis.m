function varargout = PopAnalysis(varargin)
% POPANALYSIS M-file for PopAnalysis.fig
%      POPANALYSIS, by itself, creates a new POPANALYSIS or raises the existing
%      singleton*.
%
%      H = POPANALYSIS returns the handle to a new POPANALYSIS or the handle to
%      the existing singleton*.
%
%      POPANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POPANALYSIS.M with the given input arguments.
%
%      POPANALYSIS('Property','Value',...) creates a new POPANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PopAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PopAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PopAnalysis

% Last Modified by GUIDE v2.5 08-Dec-2009 20:05:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PopAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @PopAnalysis_OutputFcn, ...
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


% --- Executes just before PopAnalysis is made visible.
function PopAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PopAnalysis (see VARARGIN)

% Choose default command line output for PopAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PopAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


global popState

popState.alignflag = 0;
popState.peakflag = 0;
popState.dFThresh = str2num(get(handles.dFThresh,'String'));
Nroi = str2num(get(handles.noROI,'string'));
popState.funcSymbol = get(handles.funcSymbol,'string');
popState.bwPopAnalysis = cell(1,Nroi);

% --- Outputs from this function are returned to the command line.
function varargout = PopAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in popTimeCourse.
function popTimeCourse_Callback(hObject, eventdata, handles)
% hObject    handle to popTimeCourse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global popState

popState.peakflag = get(handles.normPeak,'Value');
popState.alignflag = get(handles.alignPeaks,'Value');
popState.dFThresh = str2num(get(handles.dFThresh,'String'));
popState.funcSymbol = get(handles.funcSymbol,'string');

ROIcompareTime

% --- Executes on button press in popTuningCurve.
function popTuningCurve_Callback(hObject, eventdata, handles)
% hObject    handle to popTuningCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global popState

popState.peakflag = get(handles.normPeak,'Value');
popState.alignflag = get(handles.alignPeaks,'Value');
popState.dFThresh = str2num(get(handles.dFThresh,'String'));
popState.funcSymbol = get(handles.funcSymbol,'string');

ROIcompareCurve


% --- Executes on button press in normPeak.
function normPeak_Callback(hObject, eventdata, handles)
% hObject    handle to normPeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normPeak



function noROI_Callback(hObject, eventdata, handles)
% hObject    handle to noROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noROI as text
%        str2double(get(hObject,'String')) returns contents of noROI as a double


% --- Executes during object creation, after setting all properties.
function noROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in createROIs.
function createROIs_Callback(hObject, eventdata, handles)
% hObject    handle to createROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global popState

[Im] = Load2phImage(1,[1 0 0],2);

Nroi = str2num(get(handles.noROI,'string'));
popState.bwPopAnalysis = cell(1,Nroi);

colorid = 'brgky';

figure, imagesc(Im{1}), colormap gray

for i = 1:Nroi

    [popState.bwPopAnalysis{i} popState.ROIPolyx{i} popState.ROIPolyy{i}] = roipoly;
    
    hold on
    
    plot(popState.ROIPolyx{i},popState.ROIPolyy{i},colorid(i),'lineWidth',5)
 
    
end
hold off

set(handles.popTimeCourse,'enable','on')
set(handles.popTuningCurve,'enable','on')


% --- Executes on button press in alignPeaks.
function alignPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to alignPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alignPeaks



function dFThresh_Callback(hObject, eventdata, handles)
% hObject    handle to dFThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dFThresh as text
%        str2double(get(hObject,'String')) returns contents of dFThresh as a double


% --- Executes during object creation, after setting all properties.
function dFThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dFThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function funcSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to funcSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of funcSymbol as text
%        str2double(get(hObject,'String')) returns contents of funcSymbol as a double


% --- Executes during object creation, after setting all properties.
function funcSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to funcSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


