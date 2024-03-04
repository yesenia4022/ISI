pF0

%initialize the Gui

global G_handles Analyzer

% set(G_handles.epistart,'String','1000');  %Frame start in ms (to average)
% set(G_handles.epistop,'String','4000'); %Frame stop in ms (to average)
% set(G_handles.bstart,'String','-1500');  %Frame start in ms (to average)
% set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

set(G_handles.datadir,'string','F:\ISIdata\Mouse ISI\')

%%

clear ExpList

idx = 1;
ExpList{idx}.anim = 'ne8'; %nice
ExpList{idx}.Rexpt = 'u001_004'; 
ExpList{idx}.Cexpt = 'u004_001';
ExpList{idx}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{idx}.E_yx = [435 315; 380 300; 325 280];
    
idx = idx+1;
ExpList{idx}.anim = 'nf0'; %nice
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_001';
ExpList{idx}.areaNames = {'V1','AL','LM'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nf3'; %Consistent, but color map has weak dynamic range and biased to M
ExpList{idx}.Rexpt = 'u000_001'; 
ExpList{idx}.Cexpt = 'u001_006';
ExpList{idx}.areaNames = {'PM','V1','A','X','RL','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'ng4'; %awesome
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_003';
ExpList{idx}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'ng5'; %Nice color but noisy retinotopy maps. Cut off fov
ExpList{idx}.Rexpt = 'u002_001'; 
ExpList{idx}.Cexpt = 'u002_004';
ExpList{idx}.areaNames = {'PM','V1','RL','AL','LM'};
ExpList{idx}.E_yx = [435 315;380 300; 325 280];

idx = idx+1;
ExpList{idx}.anim = 'nh1'; %Good. Good color map, but S dominant.
ExpList{idx}.Rexpt = 'u000_002'; 
ExpList{idx}.Cexpt = 'u000_004';
ExpList{idx}.areaNames = {'V1','LM','X'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nh3'; %Ok.  M is really dominant for some reason
ExpList{idx}.Rexpt = 'u000_001'; 
ExpList{idx}.Cexpt = 'u000_004';
ExpList{idx}.areaNames = {'X','X','X','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nj2'; %V1 and PM are great. LM and AL are very S-dominated
ExpList{idx}.Rexpt = 'u001_001'; 
ExpList{idx}.Cexpt = 'u001_002';
ExpList{idx}.areaNames = {'PM,','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nj3'; %good.
ExpList{idx}.Rexpt = 'u000_001'; 
ExpList{idx}.Cexpt = 'u000_003';
ExpList{idx}.areaNames = {'V1','PM','RL','LM'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nj4'; 
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_002';
ExpList{idx}.areaNames = {'PM','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nk5';  %Strange color
ExpList{idx}.Rexpt = 'u001_001';
ExpList{idx}.Cexpt = 'u001_003';  
% ExpList{idx}.Rexpt = 'u000_001';
% ExpList{idx}.Cexpt = 'u000_002';  
ExpList{idx}.areaNames = {'PM','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nk6'; 
ExpList{idx}.Rexpt = 'u000_001'; 
ExpList{idx}.Cexpt = 'u000_002';
ExpList{idx}.areaNames = {'X','PM','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nl2'; %OK Mostly consistent
ExpList{idx}.Rexpt = 'u001_000'; 
ExpList{idx}.Cexpt = 'u001_001';  
ExpList{idx}.areaNames = {'V1','PM','AL','LM'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nk8'; %OK Mostly consistent
ExpList{idx}.Rexpt = 'u000_001'; 
ExpList{idx}.Cexpt = 'u000_002';  
ExpList{idx}.areaNames = {'V1','PM','X','A','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nl8'; %strange color
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_001';
ExpList{idx}.areaNames = {'PM','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nm0'; %nice
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_001';
ExpList{idx}.areaNames = {'X','V1','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nm9'; %Consistent trends, but major S bias. Kernels are a bit noisy.
ExpList{idx}.Rexpt = 'u000_002'; 
ExpList{idx}.Cexpt = 'u000_003';
ExpList{idx}.areaNames = {'V1','PM','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nn6'; 
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_001';
ExpList{idx}.areaNames = {'PM','V1','RL','LM','AL'};
ExpList{idx}.E_yx = [];

idx = idx+1;
ExpList{idx}.anim = 'nn7'; 
ExpList{idx}.Rexpt = 'u000_000'; 
ExpList{idx}.Cexpt = 'u000_001';
ExpList{idx}.areaNames = {'PM','V1','LM'};
ExpList{idx}.E_yx = [];

%%
clear DataS
expDom = 1:length(ExpList);
%expDom = 17;
for i = 1:length(expDom)
    
    anim = ExpList{expDom(i)}.anim;
    Rexpt = ExpList{expDom(i)}.Rexpt; %Retinotopy experiment
    Cexpt = ExpList{expDom(i)}.Cexpt; %Color experiment
    
    set(G_handles.loadana,'string',anim)
    set(G_handles.loadexp,'string',Cexpt)
    Gsetdirectories
    load(['c:\f0 images\' anim '_' Cexpt '_f0m'],'f0m')
    
    try
        load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end) ...
            ' from ' Cexpt(2:end)],'imstate')
    catch
        load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')
    end
    
    areaNames = ExpList{expDom(i)}.areaNames;
    E_yx = ExpList{expDom(i)}.E_yx;
    
    DataS{i} = AnalyzeBGmatrix2(imstate,E_yx,areaNames);

end

%% COMBINE ALL
pixpermm = 125;
f=100
areaOrder{1} = 'V1';
areaOrder{2} = 'PM';
areaOrder{3} = 'LM';
areaOrder{4} = 'AL';
areaOrder{5} = 'RL';
%for q = 1:length(areaOrder)
    for q = 5:5
    vertAll = [];
    cmapAll = [];
    slopesAll{q} = [];
    intAll{q} = [];
    area = areaOrder{q}
    figure
    anCount = 0;
    for i = 1:length(DataS)
        imlabel = bwlabel(DataS{i}.areaBounds,4);
        Areaid = NaN;
        for j = 1:length(DataS{i}.areaNames)
            if strcmp(area,DataS{i}.areaNames{j})
                Areaid = j;
                pixid = find(imlabel(:) == j & abs(DataS{i}.vertmap(:))<45 & abs(DataS{i}.hormap(:))<45);
                break
            end
        end        
        
        if ~isnan(Areaid) %If this mouse had this area identified
            msk = zeros(size(DataS{i}.areaBounds));
            msk(pixid) = 1;
            xdom = 0:(size(msk,2)-1)/pixpermm;
            ydom = 0:(size(msk,1)-1)/pixpermm;
            
            vdum = DataS{i}.vert_vec{Areaid};
            cdum = DataS{i}.cmap_vec{Areaid};
            
            %plot scatter of this area/animal
            %ax3 = subplot(3,length(DataS),i+length(DataS)*2);
            
            figure(98+3*q + f)
            subplot(5,5,i)
            retVScolorDist(vdum,cdum)
            id = find(~isnan(vdum.*cdum));
            [r p] = corrcoef(vdum(id),cdum(id));
            r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
            title(['r=' num2str(r) '; p=' num2str(p)]);
            ylim([0 1])
            
            %vdum = zscore(vdum);
            cdum = zscore(cdum);
            
            vertAll = [vertAll; vdum];
            cmapAll = [cmapAll; cdum];
            
            slopesAll{q} = [slopesAll{q}; DataS{i}.linfit{Areaid}(1)];
            intAll{q} = [intAll{q}; DataS{i}.linfit{Areaid}(2)];
            
            %ax1 = subplot(3,length(DataS),i);
            figure(99+3*q + f)
            subplot(5,5,i)
            imagesc(xdom,ydom,DataS{i}.vertmap, 'AlphaData',msk,[-45 45]);
            title('vertical retinotopy')
            colormap jet
            axis image
            xlabel('mm')
            
            %ax2 = subplot(3,length(DataS),i+length(DataS))
            figure(100+3*q+ f)
            subplot(5,5,i)
            imagesc(xdom,ydom,DataS{i}.cmap,  'AlphaData',msk,[0 1]);
            title('%S')
            colormap parula
            axis image
            
            anCount = anCount+1;
        end
    end
    
    figure,

    retVScolorDist(vertAll,cmapAll/4)
    id = find(~isnan(vertAll.*cmapAll));
    [r p] = corrcoef(vertAll(id),cmapAll(id));
    r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
    title([area '; ' num2str(anCount) ' mice; r=' num2str(r) '; p=' num2str(p)]);
    
    cmapAllAll{q} = cmapAll;
end

%%
figure
for i = 1:length(cmapAllAll)
    mu(i) = trimmean(cmapAllAll{i},50);
    sig(i) = std(cmapAllAll{i});
    
    subplot(5,1,i)
    hist(cmapAllAll{i},[-2:.1:2])
    xlim([-2 2])
    
    

end

figure,bar(mu)
%% COMBINE good ones
goodOnes = [1 2 3 4 6 8 9 14 16]
goodOnes = [1 2 4 6 9 14]
pixpermm = 125;
f=100
areaOrder{1} = 'V1';
areaOrder{2} = 'PM';
areaOrder{3} = 'LM';
areaOrder{4} = 'AL';
areaOrder{5} = 'RL';
%for q = 1:length(areaOrder)
    for q = 1:1
    vertAll = [];
    cmapAll = [];
    slopesAll{q} = [];
    intAll{q} = [];
    area = areaOrder{q}
    figure
    anCount = 0;
    for idum = 1:length(goodOnes)
        i = goodOnes(idum)
        imlabel = bwlabel(DataS{i}.areaBounds,4);
        Areaid = NaN;
        for j = 1:length(DataS{i}.areaNames)
            if strcmp(area,DataS{i}.areaNames{j})
                Areaid = j;
                pixid = find(imlabel(:) == j & abs(DataS{i}.vertmap(:))<45 & abs(DataS{i}.hormap(:))<45);
                break
            end
        end        
        
        if ~isnan(Areaid) %If this mouse had this area identified
            msk = zeros(size(DataS{i}.areaBounds));
            msk(pixid) = 1;
            
            xdom = 0:(size(msk,2)-1)/pixpermm;
            ydom = 0:(size(msk,1)-1)/pixpermm;
            
            vdum = DataS{i}.vert_vec{Areaid};
            cdum = DataS{i}.cmap_vec{Areaid};
            
            %plot scatter of this area/animal
           
            figure(98+3*q + f)
            ax3 = subplot(length(goodOnes),3,(idum-1)*3+1)
            retVScolorDist(vdum,cdum)
            id = find(~isnan(vdum.*cdum));
            [r p] = corrcoef(vdum(id),cdum(id));
            r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
            title(['r=' num2str(r) '; p=' num2str(p)]);
            ylim([-.1 1.1])
            ylabel('%S'),xlabel('vertical eccentricity')
            
            
            %vdum = zscore(vdum);
            cdum = zscore(cdum);
            
            vertAll = [vertAll; vdum];
            cmapAll = [cmapAll; cdum];
            
            slopesAll{q} = [slopesAll{q}; DataS{i}.linfit{Areaid}(1)];
            intAll{q} = [intAll{q}; DataS{i}.linfit{Areaid}(2)];
            
            %Make window around area
            x = 1:size(msk,2);
            y = 1:size(msk,1);
            CoMX = round(sum(x.*sum(msk)/sum(msk(:))));
            CoMY = round(sum(y'.*sum(msk,2)/sum(msk(:))));
            W = 100;
            xran = CoMX-W:CoMX+W;
            yran = CoMY-W:CoMY+W;
            xran(find(xran<1 | xran>size(msk,2))) = [];
            yran(find(yran<1 | yran>size(msk,1))) = [];
            xdom = (1:length(xran))/pixpermm;
            ydom = (1:length(yran))/pixpermm;

            %figure(99+3*q + f)
            ax1 = subplot(length(goodOnes),3,(idum-1)*3+2);
            imagesc(xdom,ydom,DataS{i}.vertmap(yran,xran), 'AlphaData',msk(yran,xran),[-45 45]);
            title('vertical retinotopy')
            colormap(ax1,jet)
            axis image
            xlabel('mm')
            colorbar

            %figure(100+3*q+ f)
            ax2 = subplot(length(goodOnes),3,(idum-1)*3+3);
            imagesc(xdom,ydom,DataS{i}.cmap(yran,xran),  'AlphaData',msk(yran,xran),[-.1 1.1]);
            title('%S')
            colormap(ax2,parula)
            axis image
            colorbar
            
            anCount = anCount+1;
        end
    end
    
    figure,

    retVScolorDist(vertAll,cmapAll/4)
    id = find(~isnan(vertAll.*cmapAll));
    [r p] = corrcoef(vertAll(id),cmapAll(id));
    r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
    title([area '; ' num2str(anCount) ' mice; r=' num2str(r) '; p=' num2str(p)]);
    
    cmapAllAll{q} = cmapAll;
end
