function Gsetdirectories

%sets data and analyzer directories based on what is in the GUI.  It also
%loads the acquisition parameters

global G_handles datadir AUE

anim = get(G_handles.loadana,'String');
expt = get(G_handles.loadexp,'String');

dir = get(G_handles.analyzedir,'String'); %partial path for analyzer file
setAnalyzerDirectory([dir anim '\' ])

loadAnalyzer(expt)

fno = 1; tno = 1;
set(G_handles.frameno,'string',num2str(fno))
set(G_handles.trialno,'string',num2str(tno))

dir = get(G_handles.datadir,'String'); %partial path for .tiff files 

datadir = [dir anim '\' expt '\'];  %Make path for .tiff files
AUE = [anim '_' expt]; %loaded unit experiment 'u000_000'

setacqinfo(tno)  %Set the global ACQinfo structure (contains trial info