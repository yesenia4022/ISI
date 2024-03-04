pF0_ISI

%initialize the Gui

global G_handles Analyzer

% set(G_handles.epistart,'String','1000');  %Frame start in ms (to average)
% set(G_handles.epistop,'String','4000'); %Frame stop in ms (to average)
% set(G_handles.bstart,'String','-1500');  %Frame start in ms (to average)
% set(G_handles.bstop,'String','0'); %Frame stop in ms (to average)
set(G_handles.basesub,'Value',1); %baseline subtraction

set(G_handles.datadir,'string','I:\ISIdata\')

%%

clear ExpList

idx = 1;
% ExpList{idx}.anim = 'ra5';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_000'; %retinotopy
% ExpList{idx}.areaNames = {'XX','V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_ra5 v2'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.Cexpt = {'u002_006','u002_002','u002_003','u002_004','u002_005'}; %place in order of highest to lowest baseline
% ExpList{idx}.E_yx = []; %ignore this
%     
% idx = idx+1;
% ExpList{idx}.anim = 'rc2';  %retinotopy
% ExpList{idx}.Rexpt = 'u002_002'; 
% ExpList{idx}.Cexpt = {'u002_003','u002_006','u002_007','u002_008','u002_012'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'XX','V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rc2'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% idx = idx+1;
% ExpList{idx}.anim = 'rb6';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_000'; 
% ExpList{idx}.Cexpt = {'u003_002','u003_005','u003_006','u003_007','u003_010'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'X','V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rb6'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rf5';  %retinotopy
% ExpList{idx}.Rexpt = 'u000_000'; 
% ExpList{idx}.Cexpt = {'u000_002','u000_005','u000_006','u000_007','u000_010'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rf5'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rf6';  %retinotopy
% ExpList{idx}.Rexpt = 'u000_000'; 
% ExpList{idx}.Cexpt = {'u000_003','u000_006','u000_007','u000_008','u000_011'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rf6'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rf4';  %retinotopy
% ExpList{idx}.Rexpt = 'u000_001'; 
% ExpList{idx}.Cexpt = {'u000_003','u000_006','u000_008','u000_009','u000_012'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rf4'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this
% 
% idx = idx+1;
% ExpList{idx}.anim = 'rf3';  %retinotopy
% ExpList{idx}.Rexpt = 'u000_002'; 
% ExpList{idx}.Cexpt = {'u000_004','u000_007','u000_008','u000_009','u000_012'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rf3'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rf7';  %retinotopy
% ExpList{idx}.Rexpt = 'u000_001'; 
% ExpList{idx}.Cexpt = {'u000_003','u000_004','u000_005','u000_008','u000_009'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rf7'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rf9';  %retinotopy
% ExpList{idx}.Rexpt = 'u000_000'; 
% ExpList{idx}.Cexpt = {'u000_001','u000_004','u000_005','u000_006','u000_009'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rf9'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rb8';  %retinotopy
% ExpList{idx}.Rexpt = 'u002_000'; 
% ExpList{idx}.Cexpt = {'u002_003','u002_005','u002_006','u002_007','u002_008'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rb8'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'nz0';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_001'; 
% ExpList{idx}.Cexpt = {'u003_003','u003_014','u003_015','u003_016','u003_022'}; %Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_nz0 v2'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = []; %ignore this

% ExpList{idx}.anim = 'rc3';  %retinotopy
% ExpList{idx}.Rexpt = 'u006_000'; 
% % ExpList{idx}.Cexpt = {'u006_010','u006_004','u006_005','u006_006','u006_009'}; %Highest to lowest baseline
% % ExpList{idx}.Cexpt = {'u006_011','u006_010','u006_004','u006_005','u006_006','u006_009'}; % V2 Highest to lowest baseline
% % ExpList{idx}.Cexpt = {'u006_014','u006_015','u006_017','u006_018','u006_019','u006_020'}; %V3 Highest to lowest baseline
% ExpList{idx}.Cexpt = {'u006_014','u006_010','u006_017','u006_018','u006_006','u006_020'}; %V3 Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rc3 v3'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rc4';  %retinotopy
% ExpList{idx}.Rexpt = 'u005_001'; 
% ExpList{idx}.Cexpt = {'u005_014','u005_003','u005_004','u005_005','u005_006','u005_009'}; %V3 Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rc3 v3'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rk3';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_020'; 
% ExpList{idx}.Cexpt = {'u003_002','u003_003','u003_006','u003_007','u003_010','u003_011'}; %V3 Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rk3'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rk4';  %retinotopy
% ExpList{idx}.Rexpt = 'u002_000'; 
% ExpList{idx}.Cexpt = {'u002_016','u002_002','u002_005','u002_014','u002_008','u002_011'}; %V3 Highest to lowest baseline
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rk4'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this


% ExpList{idx}.anim = 'rk5';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_000'; 
% ExpList{idx}.Cexpt = {'u003_015','u003_001','u003_004','u003_006','u003_007','u003_011'}; 
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rk5'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rl1';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_000'; 
% ExpList{idx}.Cexpt = {'u003_001','u003_002','u003_006','u003_007','u003_008','u003_012'}; 
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rl1'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rl0';  %retinotopy
% ExpList{idx}.Rexpt = 'u003_002'; 
% ExpList{idx}.Cexpt = {'u003_006','u003_003','u003_017','u003_018','u003_010','u003_013'}; 
% % ExpList{idx}.Cexpt = {'u003_006','u003_003','u003_007','u003_018','u003_010','u003_013'}; %B4= exp 7
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rl0'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% % ExpList{idx}.GtrxName = 'Gtrx_rl0 v7';
% ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rk6';  %retinotopy
% ExpList{idx}.Rexpt = 'u001_002'; 
% ExpList{idx}.Cexpt = {'u001_019','u001_004','u001_005','u001_006','u001_007','u001_008'}; 
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rk6'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this

ExpList{idx}.anim = 'rk7';  %retinotopy
ExpList{idx}.Rexpt = 'u002_000'; 
ExpList{idx}.Cexpt = {'u002_001','u002_002','u002_003','u002_004','u002_005','u002_006'}; 
ExpList{idx}.areaNames = {'V1'}; %area id
ExpList{idx}.GtrxName = 'Gtrx_rk7'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
ExpList{idx}.E_yx = [1]; %ignore this

% ExpList{idx}.anim = 'rl7';  %retinotopy %%%%First set of Baselines
% ExpList{idx}.Rexpt = 'u003_001'; 
% ExpList{idx}.Cexpt = {'u003_002','u003_003','u003_004','u003_005','u003_006','u003_007'}; 
% ExpList{idx}.areaNames = {'V1'}; %area id
% ExpList{idx}.GtrxName = 'Gtrx_rl7'; %file name that gets saved after running 'alignAllBaselinestoRetinotopy'
% ExpList{idx}.E_yx = [1]; %ignore this



%%
goodOnes = [1]
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
    Cexpt = ExpList{expDom(i)}.Cexpt; %Color experiments
    
    try
        load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end) ...
            ' from ' Cexpt{bL}(2:end)],'imstate')
    catch
        load(['c:\translated Kretinotopy\imstate_' anim '_' Rexpt(2:end)],'imstate')
    end
    
    
    
    for bL = 1:length(Cexpt) %Loop each baseline
        
        load(['c:\baselineTrx\' ExpList{i}.GtrxName],'trx')
        Gtrx = trx{bL};
        
        set(G_handles.loadana,'string',anim)
        set(G_handles.loadexp,'string',Cexpt{bL})
        Gsetdirectories
        load(['c:\f0 images\' anim '_' Cexpt{bL} '_f0m'],'f0m')        
        
        areaNames = ExpList{expDom(i)}.areaNames;
        E_yx = ExpList{expDom(i)}.E_yx;
        
        DataS{i,bL} = AnalyzeTheta(imstate,f0m,E_yx,Gtrx,areaNames);
        
    end
    
end


%%
for i = 1:length(expDom)
    figure
    subplot(1,length(Cexpt)+1,1)
    imagesc(DataS{i,bL}.vertmap,'AlphaData',DataS{i,bL}.areaBounds,[-50 50]), colormap jet, colorbar
    axis image
    axis off
    for bL = 1:length(Cexpt) %Loop each baseline
        subplot(1,length(Cexpt)+1,bL+1)
        imagesc(DataS{i,bL}.cmap,'AlphaData',DataS{i,bL}.areaBounds,[0 1]), colormap jet
        axis image
        axis off
    end
    colorbar
end


%%
% basedom = [.5/16 .5/8 .5/4 .5/2 .5]; %for log Baseline experiments
basedom = [.03 .1475 .265 .3825 .5 1] % for linear baseline experiments
for i = 1:length(basedom);
    legstr{i} = num2str(basedom(i));
end

clear mumat musig
for i = 1:length(expDom)

    
    for j = 1:length(ExpList{i}.areaNames)
        if strcmp(ExpList{i}.areaNames{j},'V1');
            V1id = j;
        end
    end
    
    vertdom = DataS{i,1}.ebardom{V1id};
    for v = 1:length(vertdom)
        legstr_vert{v} = num2str(round(vertdom(v)));
    end
    
    for bL = 1:6
        ebarmu = DataS{i,bL}.ebarmu{V1id};
        ebarsig = DataS{i,bL}.ebarsig{V1id};
        ebardom = DataS{i,bL}.ebardom{V1id};
        mumat(:,bL) = ebarmu;
        musig(:,bL) = ebarsig;
        
        slp(bL) = DataS{i,bL}.linfit{V1id}(2);
    end
    
    mumat = fliplr(mumat);
    musig = fliplr(musig);
    
    slp = fliplr(slp);
    figure,
    subplot(2,2,1)
    plot(ebardom,mumat,'o-'), xlabel('vertical retinotopy'), ylabel('%S')
    ylim([-.1 1])
    legend(legstr)
    subplot(2,2,2)
%     semilogx(basedom,mumat','o-'), xlabel('baseline'), ylabel('%S') %for log Baseline experiments
    plot(basedom,mumat','o-'), xlabel('baseline'), ylabel('%S') % for linear baseline experiments
    set(gca,'XTick',basedom)
    ylim([-.1 1])
    legend(legstr_vert)
    
    %Get the model
    [ret,s,base] = svd(mumat);
    ret = sign(ret(end,1))*ret; base = sign(base(end,1))*base;
    mumat_sep = ret(:,1)*base(:,1)'*s(1,1);
    
    [pRet ffit_ret varaccount] = Sigfit((ebardom),ret(:,1)',[]); 
     
    G = [median(log(basedom)) 10 (max(base(:,1))-min(base(:,1))) 0];
    [pBase ffit_base varaccount] = Sigfit(log(basedom),base(:,1)',G);
    
    %[param_base ffit_base MSE] = Expfit((basedom),1./base(:,1)',[10 .01 1]); This is the inverse of the sigmoid, w/o shift param;
    %[param_base ffit_base varaccount] = DecayFit((basedom),base(:,1)',[(.125*20) 0]);
    
    
%     basedomI = log(basedom(1)):.01:log(basedom(end));
%     ebardomI = ebardom(1):.1:ebardom(end);
%     ffit_retI = pRet(3)./(1 + exp(-(ebardomI-pRet(1))*pRet(2))) + pRet(4);
%     ffit_baseI = pBase(3)./(1 + exp(-(basedomI-pBase(1))*pBase(2))) + pBase(4);
       
    pS_model = ffit_base(:)*ffit_ret(:)'*s(1,1);
     
    subplot(2,2,3)
    plot(ebardom,pS_model','o-')   
    ylim([-.1 1])
    
    subplot(2,2,4)
%   semilogx(basedom,pS_model,'o-') %for log Baseline experiments
    plot(basedom,pS_model,'o-') % for linear baseline experiments
    set(gca,'XTick',basedom)
    xlabel('baseline'), ylabel('%S')
    ylim([-.1 1])
   
    %figure,plot(mumat_sep)
    
end





%%
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
                
                figure(100+i),
                scatter(vdum(id),cdum(id),'.k'), hold on, plot(vdum(id),lfit), hold on, scatter(vdum(id),sffit,'.r')
                
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
    
    DataS{i} = AnalyzeBGmatrix3(imstate,f0m,E_yx,areaNames);

end