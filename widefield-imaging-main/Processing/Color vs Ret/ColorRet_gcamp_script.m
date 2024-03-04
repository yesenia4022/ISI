pF0_ISI

%initialize the Gui

global G_handles Analyzer

% set(G_handles.epistart,'String','1000');  %Frame start in ms (to average)
% set(G_handles.epistop,'String','4000'); %Frame stop in ms (to average)
% set(G_handles.bstart,'String','-1500');  %Frame start in ms (to average)
% set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

set(G_handles.datadir,'string','E:\ISIdata\')

%%

clear ExpList
idx = 1;

% idx = idx;
% ExpList{idx}.anim = 'nx7'; %nice
% ExpList{idx}.Rexpt = 'u004_000'; 
% ExpList{idx}.Cexpt = 'u004_001';
% ExpList{idx}.areaNames = {'PM','V1','XX','AL','X2','LM'};
% ExpList{idx}.E_yx = [435 315; 380 300; 325 280];
%     
% idx = idx+1;
% ExpList{idx}.anim = 'nx8'; %nice
% ExpList{idx}.Rexpt = 'u004_000'; 
% ExpList{idx}.Cexpt = 'u004_002';
% ExpList{idx}.areaNames = {'PM','A','V1','LM','AL'};
% ExpList{idx}.E_yx = [];
% 
% idx = idx+1;
% ExpList{idx}.anim = 'ra2'; %ISI 
% ExpList{idx}.Rexpt = 'u000_001'; 
% ExpList{idx}.Cexpt = 'u000_002';
% ExpList{idx}.areaNames = {'V1','PM','LM'};
% ExpList{idx}.E_yx = [];
% 
% idx = idx+1;
% ExpList{idx}.anim = 'nr4'; %nice
% ExpList{idx}.Rexpt = 'u008_000'; 
% ExpList{idx}.Cexpt = 'u008_001';
% ExpList{idx}.areaNames = {'PM','A','V1','LM','AL'};
% ExpList{idx}.E_yx = [];

% %idx = idx+1;
% ExpList{idx}.anim = 'rk7';  
% ExpList{idx}.Rexpt = 'u002_000'; 
% %ExpList{idx}.Cexpt = 'u002_007'; %B5 sfreq
% % ExpList{idx}.Cexpt = 'u002_009';%B2 sfreq 
% % ExpList{idx}.Cexpt = 'u002_008';%B5 tperiod
% ExpList{idx}.Cexpt = 'u002_010';%B2 tperiod 
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% %idx = idx+1;
% ExpList{idx}.anim = 'rk6';  
% ExpList{idx}.Rexpt = 'u001_002'; 
% %ExpList{idx}.Cexpt = 'u001_009'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u001_013'; %B2 s_freq
% %ExpList{idx}.Cexpt = 'u001_010'; %B5 t_period
% ExpList{idx}.Cexpt = 'u001_014'; %B2 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% %idx = idx+1;
% ExpList{idx}.anim = 'rl0';
% ExpList{idx}.Rexpt = 'u003_002'; 
% % ExpList{idx}.Cexpt = 'u003_004'; %B5 s_freq
% % ExpList{idx}.Cexpt = 'u003_011'; %B2 s_freq
%  ExpList{idx}.Cexpt = 'u003_005'; %B5 t_period
% % ExpList{idx}.Cexpt = 'u003_012'; %B2 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% idx = idx+1;
% ExpList{idx}.anim = 'rk4';
% ExpList{idx}.Rexpt = 'u003_007'; %s_frequency
% ExpList{idx}.Cexpt = 'u003_002'; %B5- 'u003_002' B2- 'u002_009' 
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% idx = idx;
% ExpList{idx}.anim = 'rk4'; 
% ExpList{idx}.Rexpt = 'u003_007'; 
% %ExpList{idx}.Cexpt = 'u003_002'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u003_005'; %B2 s_freq
% %ExpList{idx}.Cexpt = 'u003_003'; %B5 t_period
% ExpList{idx}.Cexpt = 'u003_006'; %B2 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rm9'; 
% ExpList{idx}.Rexpt = 'u001_000'; 
% %ExpList{idx}.Cexpt = 'u001_002'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u001_008'; %B1 s_freq
% %ExpList{idx}.Cexpt = 'u001_003'; %B5 t_period
% ExpList{idx}.Cexpt = 'u001_009'; %B1 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];


% % idx = idx+1;
% ExpList{idx}.anim = 'rm9';
% ExpList{idx}.Rexpt = 'u001_000'; 
% %ExpList{idx}.Cexpt = 'u001_002'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u001_008'; %B1 s_freq
% %ExpList{idx}.Cexpt = 'u001_003'; %B5 t_period
% ExpList{idx}.Cexpt = 'u001_009'; %B1 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rm8'; 
% ExpList{idx}.Rexpt = 'u003_000'; 
% %ExpList{idx}.Cexpt = 'u003_005'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u003_009'; %B1 s_freq
% %ExpList{idx}.Cexpt = 'u003_004'; %B5 t_period
% ExpList{idx}.Cexpt = 'u003_010'; %B1 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rl5'; % unit 5 
% ExpList{idx}.Rexpt = 'u004_010'; %I can't analyze this... 
% ExpList{idx}.Cexpt = 'u005_003'; %B5- 'u005_003' B2- 'u005_011' 
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];
% ExpList{idx}.GtrxName = 'Gtrx_rl5 v1'; %B5 s_freq

% % idx = idx+1;
% ExpList{idx}.anim = 'rl5'; % unit 6
% ExpList{idx}.Rexpt = 'u006_001'; 
% %ExpList{idx}.Cexpt = 'u006_003'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u006_006'; %B1 s_freq
% %ExpList{idx}.Cexpt = 'u006_004'; %B5 t_period
% ExpList{idx}.Cexpt = 'u006_007'; %B1 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rm4'; 
% ExpList{idx}.Rexpt = 'u001_001'; 
% %ExpList{idx}.Cexpt = 'u001_003'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u001_008'; %B2 s_freq
% %%ExpList{idx}.Cexpt = 'u006_004'; %B5 t_period missing
% ExpList{idx}.Cexpt = 'u001_007'; %B2 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rm0'; 
% ExpList{idx}.Rexpt = 'u001_001'; 
% %ExpList{idx}.Cexpt = 'u001_003'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u001_006'; %B1 s_freq
% %ExpList{idx}.Cexpt = 'u001_004'; %B5 t_period
% ExpList{idx}.Cexpt = 'u001_008'; %B1 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rk5'; 
% ExpList{idx}.Rexpt = 'u003_000'; 
% ExpList{idx}.Cexpt = 'u003_002'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u003_009'; %B1 s_freq
% %ExpList{idx}.Cexpt = 'u003_003'; %B5 t_period
% %ExpList{idx}.Cexpt = 'u003_010'; %B1 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rj9'; 
% ExpList{idx}.Rexpt = 'u000_002'; 
%  ExpList{idx}.Cexpt = 'u000_004'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u000_008'; %B2 s_freq
% %ExpList{idx}.Cexpt = 'u000_005'; %B5 t_period
% %ExpList{idx}.Cexpt = 'u000_09'; %B2 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rn4'; 
% ExpList{idx}.Rexpt = 'u004_001'; 
% %ExpList{idx}.Cexpt = 'u004_004'; %B5 s_freq
% ExpList{idx}.Cexpt = 'u004_007'; %B2 s_freq
% %ExpList{idx}.Cexpt = 'u000_000'; %B5 t_period
% %ExpList{idx}.Cexpt = 'u000_009'; %B2 t_period
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rm1'; 
% ExpList{idx}.Rexpt = 'u001_000'; 
% %ExpList{idx}.Cexpt = 'u001_002'; %B5 s_freq
% ExpList{idx}.Cexpt = 'u001_004'; %B1 s_freq
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'rl7'; 
% ExpList{idx}.Rexpt = 'u005_000'; 
% %ExpList{idx}.Cexpt = 'u005_002'; %B5 s_freq
% %ExpList{idx}.Cexpt = 'u001_004'; %B1 s_freq
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% % idx = idx+1;
% ExpList{idx}.anim = 'se3'; 
% ExpList{idx}.Rexpt = 'u001_002'; 
% % ExpList{idx}.Cexpt = 'u001_004'; %B5 [0 90] s_freq
% % ExpList{idx}.Cexpt = 'u001_005'; %B5 [45 135] s_freq
% ExpList{idx}.areaNames = {'V1'};
% ExpList{idx}.E_yx = [];

% idx = idx+1;
ExpList{idx}.anim = 'se4'; 
ExpList{idx}.Rexpt = 'u001_001'; 
% ExpList{idx}.Cexpt = 'u001_003'; %B5 [0 90] s_freq
ExpList{idx}.Cexpt = 'u001_004'; %B5 [45 135] s_freq
ExpList{idx}.areaNames = {'V1'};
ExpList{idx}.E_yx = [];


%%
goodOnes = [1]
pixpermm = 65;

clear DataS

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
    
    figure,imagesc(imstate.fmaps{1})
    bw = roipoly;
    imstate.areaBounds = bw;
    
    areaNames = ExpList{expDom(i)}.areaNames;
    E_yx = ExpList{expDom(i)}.E_yx;
    
%   DataS{i} = AnalyzeColorTFmatrix2(imstate,E_yx,areaNames);
%     DataS{i} = AnalyzeColorSFmatrix(imstate,E_yx,areaNames);
    DataS{i} = AnalyzeColor_090_SFmatrix(imstate,E_yx,areaNames); % for SF done with theta broken up into [0 90] and [45 135]


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