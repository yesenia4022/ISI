 function varargout = OverlayGuide(varargin)
% OVERLAYGUIDE M-file for OverlayGuide.fig
%      OVERLAYGUIDE, by itself, creates a new OVERLAYGUIDE or raises the existing
%      singleton*.
%
%      H = OVERLAYGUIDE returns the handle to a new OVERLAYGUIDE or the handle to
%      the existing singleton*.
%
%      OVERLAYGUIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OVERLAYGUIDE.M with the given input arguments.
%
%      OVERLAYGUIDE('Property','Value',...) creates a new OVERLAYGUIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OverlayGuide_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OverlayGuide_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OverlayGuide

% Last Modified by GUIDE v2.5 05-Jun-2016 12:08:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OverlayGuide_OpeningFcn, ...
                   'gui_OutputFcn',  @OverlayGuide_OutputFcn, ...
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


% --- Executes just before OverlayGuide is made visible.
function OverlayGuide_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OverlayGuide (see VARARGIN)

% Choose default command line output for OverlayGuide
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OverlayGuide wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global PF1handles imstate Gtimer

PF1handles = varargin{1};

reset_imstate

%construct a timer
Gtimer = timer;
set(Gtimer,'Period',0.5,'BusyMode','drop','ExecutionMode',...
    'fixedSpacing','TimerFcn',@GrabContOverlay)

if ~isfield(imstate,'imanat')
    set(handles.intensitySlider,'enable','off')
    set(handles.flipUpDown,'enable','off')
    set(handles.flipLeftRight,'enable','off')
    set(handles.rotateImage,'enable','off')
end

% --- Outputs from this function are returned to the command line.
function varargout = OverlayGuide_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in flipLeftRight.
function flipLeftRight_Callback(hObject, eventdata, handles)
% hObject    handle to flipLeftRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resetImage('LR')

if ishandle(86)
    funcAnatomyPlot
end

% --- Executes on button press in flipUpDown.
function flipUpDown_Callback(hObject, eventdata, handles)
% hObject    handle to flipUpDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resetImage('UD')

if ishandle(86)
    funcAnatomyPlot
end

% --- Executes on button press in rotateImage.
function rotateImage_Callback(hObject, eventdata, handles)
% hObject    handle to rotateImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resetImage('rotate')

if ishandle(86)
    funcAnatomyPlot
end


% --- Executes on slider movement.
function intensitySlider_Callback(hObject, eventdata, handles)
% hObject    handle to intensitySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global imstate

imstate.intRatio = get(handles.intensitySlider,'value');

resetImage('')

if ishandle(86)
    funcAnatomyPlot
end

% --- Executes during object creation, after setting all properties.
function intensitySlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to intensitySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in loadImage.
function loadImage_Callback(hObject, eventdata, handles)
% hObject    handle to loadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global imstate Analyzer

reset_imstate


animal = get(findobj('Tag','loadana'),'String');
%unit   = get(findobj('Tag','unittxt'),'String');
%expt   = get(findobj('Tag','expttxt'),'String');
datadir= get(findobj('Tag','datadir'),'String');
%tag    = get(findobj('Tag','tagtxt'),'String');

dd = [datadir '\' lower(animal) '\grabs\'];

[filename, pathname] = uigetfile([dd '\*.mat'], 'Select a grab');

S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S

if isfield(S,'grab')
    im = S.grab.img;
else
    im = S.im;
end

%Cropping AVT works differently.  It is done at the camera, not in matlab.
%i.e. only crop the recent grab if you were using the Dalsa and set an ROI.
%Otherwise, leave it alone.
if isempty(Analyzer.ACQ.ROIcrop) | strcmp(Analyzer.ACQ.camera,'AVT')
    yran = 1:size(im,1);
    xran = 1:size(im,2);
else
    yran = Analyzer.ACQ.ROIcrop(2):(Analyzer.ACQ.ROIcrop(4)+Analyzer.ACQ.ROIcrop(2)-1);
    xran = Analyzer.ACQ.ROIcrop(1):(Analyzer.ACQ.ROIcrop(3)+Analyzer.ACQ.ROIcrop(1)-1);
end

%%%%% Panda hack %%%%%%%
xran = xran - (xran(1,1)-1);
yran = yran - (yran(1,1)-1);

im = im(yran,xran);
imstate.imanat = double(im);

imstate.intRatio = .5;  %func/anatomy weight of image (scalar)
set(handles.intensitySlider,'value',imstate.intRatio)

set(handles.intensitySlider,'enable','on')
set(handles.flipUpDown,'enable','on')
set(handles.flipLeftRight,'enable','on')
set(handles.rotateImage,'enable','on')

funcAnatomyPlot

% --- Executes on button press in plotStatic.
function plotStatic_Callback(hObject, eventdata, handles)
% hObject    handle to plotStatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

funcAnatomyPlot

% --- Executes on button press in grabContOverlay.
function grabContOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to grabContOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Gtimer h

start(Gtimer)

% --- Executes on button press in grabImage.
function grabImage_Callback(hObject, eventdata, handles)
% hObject    handle to grabImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% reset_imstate
% 
% global imstate h IMGSIZE ROIcrop imagerhandles;
% 
% h.mildig.set('GrabFrameEndEvent',0,'GrabEndEvent',...
%             0,'GrabStartEvent',0);
% 
% img = zeros(ROIcrop(3),ROIcrop(4),'uint16');
% zz  = zeros(ROIcrop(3),ROIcrop(4),'uint16');
% img = h.milimg.Get(zz,IMGSIZE^2,-1,ROIcrop(1),ROIcrop(2),ROIcrop(3),ROIcrop(4));
% 
% grab.img = img;       %% image
% grab.clock = clock;   %% time stamp
% figure(10);
% imagesc(grab.img'),axis off, colormap gray; truesize
% r = questdlg('Do you want to save this grab?','Single Grab','Yes','No','Yes');
% if(strcmp(r,'Yes'))
%     
%     grab.comment = inputdlg('Please enter description:','Image Grab',1,{'No description'},'on');
% 
%     animal = get(findobj('Tag','animaltxt'),'String');
%     unit   = get(findobj('Tag','unittxt'),'String');
%     expt   = get(findobj('Tag','expttxt'),'String');
%     datadir= get(findobj('Tag','datatxt'),'String');
%     tag    = get(findobj('Tag','tagtxt'),'String');
% 
%     dd = [datadir '\' lower(animal) '\grabs\'];
%     if(~exist(dd))
%         mkdir(dd);
%     end
%     fname = [dd 'grab_' lower(get(imagerhandles.animaltxt,'String')) '_' ...
%         get(imagerhandles.unittxt,'String') '_' ...
%         get(imagerhandles.expttxt,'String') '_' ...
%         datestr(now)];
%     fname = strrep(fname,' ','_');
%     fname = strrep(fname,':','_');
%     fname = strrep(fname,'-','_');
%     fname = [fname '.mat'];
%     fname(2) = ':';
%     
%     uisave('grab',fname);
% end
% delete(10);
% 
% 
% imstate.imanat = img;
% 
% set(handles.intensitySlider,'enable','on')
% set(handles.flipUpDown,'enable','on')
% set(handles.flipLeftRight,'enable','on')
% set(handles.rotateImage,'enable','on')


% --- Executes on button press in horRet.
function horRet_Callback(hObject, eventdata, handles)
% hObject    handle to horRet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of horRet

global imstate

flag = get(handles.vertRet,'value');

if ~flag  %if this button was already pressed don't do anything
    set(handles.horRet,'value',1);
else
    imstate.imfunc = imstate.fmaps{1};;
    set(handles.vertRet,'value',0);
    funcAnatomyPlot
end

% --- Executes on button press in vertRet.
function vertRet_Callback(hObject, eventdata, handles)
% hObject    handle to vertRet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vertRet

global imstate

flag = get(handles.horRet,'value');

if ~flag  %if this button was already pressed don't do anything
    set(handles.vertRet,'value',1);
else
    imstate.imfunc = imstate.fmaps{2};
    set(handles.horRet,'value',0);
    funcAnatomyPlot
end

% --- Executes on button press in stopGrabbing.
function stopGrabbing_Callback(hObject, eventdata, handles)
% hObject    handle to stopGrabbing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Gtimer

stop(Gtimer)


% --- Executes on button press in sigMag.
function sigMag_Callback(hObject, eventdata, handles)
% hObject    handle to sigMag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sigMag

global imstate kmap_vert f1m

val = get(handles.sigMag,'value');

if val
    imstate.mag = abs(f1m{1}) + abs(f1m{2});
else
    imstate.mag = imstate.bw;
end

funcAnatomyPlot


% --- Executes on button press in alignMapsToVasc.
function alignMapsToVasc_Callback(hObject, eventdata, handles)
% hObject    handle to alignMapsToVasc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%

global imstate kmap_hor kmap_vert f1m bw areaBounds

%Note that if you flipped or rotated this will reset stuff
imstate.fmaps{1} = kmap_hor;
imstate.fmaps{2} = kmap_vert;
imstate.sigMag = abs(f1m{1}) + abs(f1m{2});
imstate.bw = bw;
imstate.imfunc = kmap_hor;
imstate.mag = bw;
%imstate.imanat = ones(size(bw));
imstate.imanatFunc = getTrialFrame(1,1);
imstate.areaBounds = areaBounds;

r = questdlg('Do you want select new input and output images for transformation?','Yes','No');
if(strcmp(r,'Yes'))
    
    dd = 'c:\';
    
    [filename, pathname] = uigetfile([dd '/*.mat'], 'Select "Input" image');
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    if isfield(S,'grab')
        im_input = double(S.grab.img);
    else
        im_input = double(S.im);
    end
    
    [filename, pathname] = uigetfile([dd '/*.mat'], 'Select "Output" image');
    S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
    if isfield(S,'grab')
        im_output = double(S.grab.img);
    else
        im_output = double(S.im);
    end    
    
else
    
    im_input = imstate.imanatFunc; %a single frame from the functional experiment
    im_output = imstate.imanat; %The loaded image from the Overlay Guide
    
end

%Im = imstate.imanatFunc; 
im_input = im_input-min(im_input(:));
im_input = im_input/max(im_input(:));

%Imvasc = imstate.imanat; %The loaded image from the Overlay Guide
im_output = im_output-min(im_output(:));
im_output = im_output/max(im_output(:));

[movingPoints,fixedPoints] = cpselect(im_input,im_output,'Wait',true);

trx = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
Ref_Im = imref2d(size(im_input));

%%
idNaN = find(isnan(imstate.fmaps{1}));
imstate.fmaps{1}(idNaN) = 0;
idNaN = find(isnan(imstate.fmaps{2}));
imstate.fmaps{2}(idNaN) = 0;
idNaN = find(isnan(imstate.sigMag));
imstate.sigMag(idNaN) = 0;
idNaN = find(isnan(imstate.imfunc));
imstate.imfunc(idNaN) = 0;
idNaN = find(isnan(imstate.imanatFunc));
imstate.imanatFunc(idNaN) = 0;
idNaN = find(isnan(imstate.imanatFunc));
imstate.imanatFunc(idNaN) = 0;
idNaN = find(isnan(imstate.mag));
imstate.mag(idNaN) = 0;
%%

imstate.fmaps{1} = imwarp(imstate.fmaps{1},trx,'OutputView',Ref_Im);
imstate.fmaps{2} = imwarp(imstate.fmaps{2},trx,'OutputView',Ref_Im);
imstate.sigMag = imwarp(imstate.sigMag,trx,'OutputView',Ref_Im);
imstate.bw = imwarp(imstate.bw,trx,'OutputView',Ref_Im);
imstate.imfunc = imwarp(imstate.imfunc,trx,'OutputView',Ref_Im);
imstate.imanatFunc = imwarp(imstate.imanatFunc,trx,'OutputView',Ref_Im);
imstate.mag = imwarp(imstate.mag,trx,'OutputView',Ref_Im);
%imstate.imanat = imwarp(imstate.imanat,trx,'OutputView',Ref_Im);

if ~isempty(imstate.areaBounds)
    imstate.areaBounds = floor(imwarp(imstate.areaBounds,trx,'OutputView',Ref_Im));
    SE = strel('disk',1);
    imstate.areaBounds = imopen(imstate.areaBounds,SE); %Trx tends to fuse areas
end

Im_registered = imwarp(im_input,trx,'OutputView',Ref_Im);

figure,
subplot(1,3,1), imagesc(im_input), title('unregistered input image'), axis image
for i =1:length(fixedPoints(:,1))
    hold on,
    plot(fixedPoints(i,1),fixedPoints(i,2),'xr')
end
    
    
subplot(1,3,2), imagesc(Im_registered), title('registered input image'), axis image
for i =1:length(fixedPoints(:,1))
    hold on,
    plot(fixedPoints(i,1),fixedPoints(i,2),'xr')
end

subplot(1,3,3), imagesc(im_output), title('output image'), axis image
for i =1:length(fixedPoints(:,1))
    hold on,
    plot(fixedPoints(i,1),fixedPoints(i,2),'xr')
end
colormap gray


funcAnatomyPlot




% --- Executes on button press in save_imstate.
function save_imstate_Callback(hObject, eventdata, handles)
% hObject    handle to save_imstate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global imstate Analyzer

UE = [Analyzer.M.unit '_' Analyzer.M.expt];
path = 'C:\';
filename = strcat(path,'imstate_',Analyzer.M.anim,'_',UE);
uisave('imstate',filename)
