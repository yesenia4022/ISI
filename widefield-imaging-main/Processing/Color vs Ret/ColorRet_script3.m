pF0

%initialize the Gui

global G_handles Analyzer

% set(G_handles.epistart,'String','1000');  %Frame start in ms (to average)
% set(G_handles.epistop,'String','4000'); %Frame stop in ms (to average)
% set(G_handles.bstart,'String','-1500');  %Frame start in ms (to average)
% set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

set(G_handles.datadir,'string','G:\ISIdata\Mouse ISI\')

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
goodOnes = [1 2 3 4 6 8 9 11 13 14 16 17]
%goodOnes = [1 2 4 6 9 14]
pixpermm = 125;

%goodOnes = [8]
clear DataS
%goodOnes = [1 2 3 4 6 8 9 11 13 14 16 17]

expDom = goodOnes

clear DataS
%expDom = 1:length(ExpList);
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

%% COMBINE ALL (This generates fig 4 of the paper)
[Sper] = getRetinaGradient(-45:45)
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
    sigmoidSlopes{q} = [];
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
            [mapmudum mapsigdum vertDomaindum] = retVScolorDist(vdum,cdum);
            id = find(~isnan(vdum.*cdum));
            [r p] = corrcoef(vdum(id),cdum(id));
            r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
            title(['r=' num2str(r) '; p=' num2str(p)]);
            ylim([0 1])
            
            if q == 1
                Emax = 45;
                sigdom = -Emax:Emax;                
                param = DataS{i}.sigfit{Areaid};
                d = sigdom-param(1);
                sigffit = param(3)*(1./(1 + exp(-(d)*param(2))) - .5) + param(4);
                slopeAtorigin = round(param(3)*param(2)/4*1000)/10; %perc/deg
                sigmoidSlopes{q} = [sigmoidSlopes{q} slopeAtorigin];
                hold on, plot(sigdom,sigffit,'b')
                base = DataS{i}.linfit{Areaid}(2);
                sigmax(i) = max(sigffit);
                sigmin(i) = min(sigffit);
                
                mapmu{i} = mapmudum;
                mapsig{i} = mapsigdum;
                vertDomain{i} = vertDomaindum;
                
                id = find(vdum>-45 & vdum<45 & cdum>-.1 & cdum<1.5);
                lfit = vdum(id)*DataS{i}.linfit{Areaid}(1) + DataS{i}.linfit{Areaid}(2);
                sffit = param(3)*(1./(1 + exp(-(vdum(id)-param(1))*param(2))) - .5) + param(4);
                sigvaracc(i) = (var(cdum(id))-(var(cdum(id)-sffit)))/var(cdum(id));
                linvaracc(i) = (var(cdum(id))-(var(cdum(id)-lfit)))/var(cdum(id));
                relvar(i) = var(cdum(id)-sffit)/var(cdum(id)-lfit);
                
                figure(100+i),scatter(vdum(id),cdum(id),'.k'), hold on, plot(vdum(id),lfit), hold on, scatter(vdum(id),sffit,'.r')
                
            end
            
            
            %vdum = zscore(vdum);
            %cdum = zscore(cdum);
            vdum = vdum-median(vdum);
            cdum = cdum-median(cdum);
            %cdum = cdum-DataS{i}.linfit{Areaid}(2);
           % vdum = vdum- (-DataS{i}.linfit{Areaid}(2)/DataS{i}.linfit{Areaid}(1));
            
            vertAll = [vertAll; vdum(id)];
            cmapAll = [cmapAll; cdum(id)];
            
            slopesAll{q} = [slopesAll{q}; DataS{i}.linfit{Areaid}(1)*100]; % perc/deg
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
    
    
    H = vertAll;
    slp = inv(H'*H)*H'*cmapAll;
    slp = round(slp*10000)/100; %perc/deg
    
    figure(222),    
    subplot(5,1,q)
    retVScolorDist(vertAll,cmapAll)
    id = find(~isnan(vertAll.*cmapAll));
    [r p] = corrcoef(vertAll(id),cmapAll(id));
    r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
    title([area '; ' num2str(anCount) ' mice; slp=' num2str(slp) '; p=' num2str(p)]);
    ylim([-.5 .5])
    set(gca,'XTick',[-50 0 50])
    xlim([-50 50])
    set(gca,'YTick',[-.50 0 .50])
    cmapAllAll{q} = cmapAll;
    if q == 1
        Emax = 45;
        [param cmapHat varaccount] = Sigfit(vertAll(id),cmapAll(id));
        sigdom = -Emax:Emax;
        d = sigdom-param(1);
        sigffit = param(3)*(1./(1 + exp(-d*param(2))) - .5) + param(4);
        slopeAtorigin = round(param(3)*param(2)/4*1000)/10; %perc/deg
        hold on, plot(sigdom,sigffit,'b')
        hold on, plot(-45:45,Sper/100-.5)
    end
    
    
end
%%
for i = 1:length(mapmu)

    
   first(i) =mapmu{i}(1);
   last(i) =mapmu{i}(end);
   smallest(i) = min(mapmu{i});
   biggest(i) = max(mapmu{i});
    
end
%%
%for i = 1:length(slopesAll)
for i = 1:5
    
    dum = slopesAll{i};
    id = find(dum>10*mean(dum) | dum<0);
    dum(id) = NaN;
    
    nanstd(dum)
    
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
goodOnes = [1 2 3 4 6 8 9 11 13 14 16 17]
goodOnes = [1 2 4 6 9 14]
pixpermm = 125;

%goodOnes = [4]
clear DataS
%goodOnes = [1 2 3 4 6 8 9 11 13 14 16 17]

expDom = goodOnes
%expDom = 17;1
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
    
    DataS{i} = AnalyzeBGmatrix3(imstate,E_yx,areaNames);

end