pF0

%initialize the Gui

global G_handles Analyzer

% set(G_handles.epistart,'String','1000');  %Frame start in ms (to average)
% set(G_handles.epistop,'String','4000'); %Frame stop in ms (to average)
% set(G_handles.bstart,'String','-1500');  %Frame start in ms (to average)
% set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

set(G_handles.datadir,'string','g:\ISIdata\Mouse ISI\')
%%
anim = 'ne8';  
Rexpt = 'u001_004'; %Retinotopy experiment
Cexpt = 'u004_001'; %Color experiment

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end) ...
    ' from ' Cexpt(2:end)],'imstate')

areaNames = {'PM','V1','RL','LM','AL'};
E_yx = [435 315; 380 300; 325 280];

DataS{1} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{1}.anim = 'ne8'; 
ExpList{1}.Rexpt = 'u001_004'; 
ExpList{1}.Cexpt = 'u004_001';
ExpList{1}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{1}.E_yx = [435 315; 380 300; 325 280];
    
%%
anim = 'nf0';
Rexpt = 'u000_000'; %Retinotopy experiment
Cexpt = 'u000_001'; %Color experiment

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end) ...
    ' from ' Cexpt(2:end)],'imstate')

areaNameVec = {'V1','X','X','X','X'};
E_yx = [];
DataS{2} = AnalyzeBGmatrix(imstate,E_yx,areaNameVec);

ExpList{2}.anim = 'nf0'; 
ExpList{2}.Rexpt = 'u000_000'; 
ExpList{2}.Cexpt = 'u000_001';
ExpList{2}.areaNames = {'V1','X','X','X','X'};
ExpList{2}.E_yx = [];

%%
anim = 'ng4'; %nice HVA maps
Rexpt = 'u000_000'; %Retinotopy experiment
Cexpt = 'u000_003'; %Color experiment

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')
areaNames = {'PM','V1','RL','LM','AL'};
E_yx = [];
DataS{3} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{3}.anim = 'ng4'; 
ExpList{3}.Rexpt = 'u000_000'; 
ExpList{3}.Cexpt = 'u000_003';
ExpList{3}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{3}.E_yx = [];
%%
% anim = 'ng5';
% Rexpt = 'u000_000'; %Retinotopy experiment
% Cexpt = 'u001_003'; %Color experiment
% 
% load(['C:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end) ...
%     ' from ' Cexpt(2:end)],'imstate')


anim = 'ng5';  %Color maps look good, but retinotopy sucks
Rexpt = 'u002_001'; %Retinotopy experiment  (vertical is weak)
Cexpt = 'u002_004'; %Color experiment

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['C:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')

load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

% load(['C:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end) ...
%     ' from ' Cexpt(2:end)],'imstate')

areaNames = {'PM','V1','RL','AL','LM'};
E_yx = [435 315;380 300; 325 280];
DataS{4} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{4}.anim = 'ng5'; 
ExpList{4}.Rexpt = 'u002_001'; 
ExpList{4}.Cexpt = 'u002_004';
ExpList{4}.areaNames = {'PM','V1','RL','AL','LM'};
ExpList{4}.E_yx = [435 315;380 300; 325 280];

%%
anim = 'nj3';
Rexpt = 'u000_001'; %Retinotopy experiment
Cexpt = 'u000_003'; %Color experiment

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')

load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')
E_yx = [];
areaNames = {'V1','PM','RL','LM'};

DataS{5} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{5}.anim = 'nj3'; 
ExpList{5}.Rexpt = 'u000_001'; 
ExpList{5}.Cexpt = 'u000_003';
ExpList{5}.areaNames = {'V1','PM','RL','LM'};
ExpList{5}.E_yx = [];
%%
anim = 'nj4'; %Something strange here.  The M and S weights are incosistent with other experiments
Rexpt = 'u000_000'; %Retinotopy experiment
Cexpt = 'u000_002'; %Color experiment


set(G_handles.datadir,'string','g:\ISIdata\Mouse ISI\')
set(G_handles.analyzedir,'string','C:\AnalyzerFiles\')

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

areaNames = {'PM','V1','LM','AL'};
E_yx = [];
DataS{6} = AnalyzeBGmatrix(imstate,E_yx,areaNames);


ExpList{6}.anim = 'nj4'; 
ExpList{6}.Rexpt = 'u000_000'; 
ExpList{6}.Cexpt = 'u000_002';
ExpList{6}.areaNames = {'PM','V1','LM','AL'};
ExpList{6}.E_yx = [];
%%
anim = 'nk5'; %Something strange here.  The M and S weights are incosistent with other experiments
Rexpt = 'u001_001'; %Retinotopy experiment
Cexpt = 'u001_003'; %Color experiment


set(G_handles.datadir,'string','f:\ISIdata\Mouse ISI\')
set(G_handles.analyzedir,'string','C:\AnalyzerFiles\')

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

areaNames = {'PM','V1','LM','AL'};
E_yx = [];
DataS{7} = AnalyzeBGmatrix(imstate,E_yx,areaNames);


ExpList{7}.anim = 'nk5'; 
ExpList{7}.Rexpt = 'u001_001'; 
ExpList{7}.Cexpt = 'u001_003';
ExpList{7}.areaNames = {'PM','V1','LM','AL'};
ExpList{7}.E_yx = [];
%%
anim = 'nk6'; 
Rexpt = 'u000_001'; %Retinotopy experiment
Cexpt = 'u000_002'; %Color experiment


set(G_handles.datadir,'string','f:\ISIdata\Mouse ISI\')
set(G_handles.analyzedir,'string','C:\AnalyzerFiles\')

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

areaNames = {'X','PM','V1','LM','AL'};
E_yx = [];
DataS{8} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{8}.anim = 'nk6'; 
ExpList{8}.Rexpt = 'u000_001'; 
ExpList{8}.Cexpt = 'u000_002';
ExpList{8}.areaNames = {'X','PM','V1','LM','AL'};
ExpList{8}.E_yx = [];
%%
anim = 'nl8'; %strange color
Rexpt = 'u000_000'; %Retinotopy experiment
Cexpt = 'u000_001'; %Color experiment


set(G_handles.datadir,'string','f:\ISIdata\Mouse ISI\')
set(G_handles.analyzedir,'string','C:\AnalyzerFiles\')

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

areaNames = {'PM','V1','LM','AL'};
E_yx = [];
DataS{9} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{9}.anim = 'nl8'; %strange color
ExpList{9}.Rexpt = 'u000_000'; 
ExpList{9}.Cexpt = 'u000_001';
ExpList{9}.areaNames = {'PM','V1','LM','AL'};
ExpList{9}.E_yx = [];

%%
anim = 'nm0'; %Nice
Rexpt = 'u000_000'; %Retinotopy experiment
Cexpt = 'u000_001'; %Color experiment


set(G_handles.datadir,'string','f:\ISIdata\Mouse ISI\')
set(G_handles.analyzedir,'string','C:\AnalyzerFiles\')

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

areaNames = {'X','V1','LM','AL'};
E_yx = [];
DataS{10} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{10}.anim = 'nm0'; nice
ExpList{10}.Rexpt = 'u000_000'; 
ExpList{10}.Cexpt = 'u000_001';
ExpList{10}.areaNames = {'X','V1','LM','AL'};
ExpList{10}.E_yx = [];

%%
anim = 'nl2'; %Color map is totally messed up
Rexpt = 'u001_000'; %Retinotopy experiment
Cexpt = 'u001_001'; %Color experiment


set(G_handles.datadir,'string','f:\ISIdata\Mouse ISI\')
set(G_handles.analyzedir,'string','C:\AnalyzerFiles\')

set(G_handles.loadana,'string',anim)
set(G_handles.loadexp,'string',Cexpt)
Gsetdirectories
load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')


load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')

areaNames = {'PM','V1','LM','AL'};
E_yx = [];
DataS{11} = AnalyzeBGmatrix(imstate,E_yx,areaNames);

ExpList{11}.anim = 'nl2'; 
ExpList{11}.Rexpt = 'u001_000'; 
ExpList{11}.Cexpt = 'u001_001';  %Color map is totally messed up
ExpList{11}.areaNames = {'PM','V1','LM','AL'};
ExpList{11}.E_yx = [];

%%

ExpList{1}.anim = 'ne8'; ExpList{1}.Rexpt = 'u001_004'; ExpList{1}.Cexpt = 'u004_001';
ExpList{1}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{1}.E_yx = [435 315; 380 300; 325 280];

ExpList{2}.anim = 'ne8'; ExpList{1}.Rexpt = 'u001_004'; ExpList{1}.Cexpt = 'u004_001';
ExpList{2}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{2}.E_yx = [435 315; 380 300; 325 280];


%% COMBINE ALL
pixpermm = 125;
vertAll = [];
cmapAll = [];
area = 'PM';
figure
for i = 1:length(DataS)   
    imlabel = bwlabel(DataS{i}.areaBounds,4);
    for j = 1:length(DataS{i}.areaNames)        
        if strcmp(area,DataS{i}.areaNames{j})
            Areaid = j;
            pixid = find(imlabel(:) == j & abs(DataS{i}.vertmap(:))<45 & abs(DataS{i}.hormap(:))<45);
            break
        end        
    end
    
    msk = zeros(size(DataS{i}.areaBounds));
    msk(pixid) = 1;
    xdom = 0:(size(msk,2)-1)/pixpermm;
    ydom = 0:(size(msk,1)-1)/pixpermm;
    
    vertAll = [vertAll; DataS{i}.vert_vec{Areaid}];
    cmapAll = [cmapAll; DataS{i}.cmap_vec{Areaid}];
    
    ax3 = subplot(3,length(DataS),i+length(DataS)*2)
    retVScolorDist(DataS{i}.vert_vec{Areaid}-mean(DataS{i}.vert_vec{Areaid}),DataS{i}.cmap_vec{Areaid})
    
    ax1 = subplot(3,length(DataS),i)    
    imagesc(xdom,ydom,DataS{i}.vertmap, 'AlphaData',msk,[-45 45]); 
    title('vertical retinotopy')
    colormap(ax1,jet)
    axis image
    xlabel('mm')
    
    ax2 = subplot(3,length(DataS),i+length(DataS))    
    imagesc(xdom,ydom,DataS{i}.cmap, 'AlphaData',msk);
    title('S-M')
    colormap(ax2,parula)
    axis image
    
end

figure,retVScolorDist(vertAll,cmapAll)
