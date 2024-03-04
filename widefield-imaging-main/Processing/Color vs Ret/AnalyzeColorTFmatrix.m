function DataS = AnalyzeColorTFmatrix(imstate,E_yx,varargin)

global f0m Analyzer pixpermm

stim_framerate = 60;

odom = getdomain('ori');
cdom = getdomain('theta');
tdom = getdomain('t_period');

imstate.areaBounds(find(imstate.areaBounds<1)) = 0;
imstate.areaBounds(find(~imstate.bw)) = 0;

%Smooth the boundary edges
sig = 1;
hh = fspecial('gaussian',size(imstate.areaBounds),sig);
imstate.areaBounds = ifft2(fft2(imstate.areaBounds).*abs(fft2(hh)));
imstate.areaBounds = (sign(imstate.areaBounds-.5)+1)/2;

DataS.areaNames = varargin{1};

pixpermm = 125;

%% Build cell matrix of responses to each color and t-period combination

imMat = cell(length(tdom),length(cdom));
for cond = 1:length(f0m)-1  
    
    theta = getVal('theta',cond);
    cid = find(cdom == theta);
    
    t_period = getVal('t_period',cond);
    tid = find(tdom == t_period);
    
    if isempty(imMat{tid,cid})
        imMat{tid,cid} = 0;
    end
    
    imMat{tid,cid} = imMat{tid,cid}+f0m{cond}/length(odom);
    
    col(tid,cid) = theta;
    tp(tid,cid) = t_period;
    
end

%spatial Filter
sig = 2;
hh = fspecial('gaussian',size(imMat{1,1}),sig);


%% Get color map and filter
colorF0 = zeros(size(imMat{1,1},1),size(imMat{1,1},2),length(cdom));
for t = 1:length(tdom)
    for c = 1:length(cdom)
        colorF0(:,:,c) = colorF0(:,:,c)+imMat{t,c}/length(tdom);
    end
end
for c = 1:length(cdom)
    colorF0(:,:,c) = ifft2(abs(fft2(hh)).*fft2(colorF0(:,:,c))); 
end

colorMap = zeros(size(imMat{1,1}));
for k = 1:length(cdom)
    colorMap = colorMap + colorF0(:,:,k)*exp(1i*2*cdom(k)*pi/180);    %Linear combination
end

colormapmag = abs(colorMap);
colorMap = angle(colorMap)/2*180/pi;
id = find(colorMap<0);
colorMap(id) = colorMap(id)+180;

%% Get a single tf tuning curve for each color, and plot

id = find(imstate.areaBounds);

tempAll = cell(1,length(cdom));
for c = 1:length(cdom)
    tempAll{c} = zeros(size(imMat{1,1},1),size(imMat{1,1},2),length(tdom));
    for t = 1:length(tdom)
         dum = imMat{t,c};
         dum = ifft2(abs(fft2(hh)).*fft2(dum));
         tempAll{c}(:,:,t) = dum;
         
         tftc(c,t) = mean(dum(id));
         
    end
    cleg{c} = num2str(cdom(c));
end

figure,semilogx(60./tdom,tftc','-o')
xlabel('t freq')
set(gca,'XTick',sort(60./tdom))

legend(cleg) 

%% Get tfreq map and filter

tempF0 = zeros(size(imMat{1,1},1),size(imMat{1,1},2),length(tdom));
for t = 1:length(tdom)
    for c = 1:length(cdom)
        tempF0(:,:,t) = tempF0(:,:,t)+imMat{t,c}/length(cdom);
    end
end
for t = 1:length(tdom)
    tempF0(:,:,t) = ifft2(abs(fft2(hh)).*fft2(tempF0(:,:,t))); 
end



% mi = min(tempF0,[],3);
% for t = 1:length(tdom)
%     tempF0(:,:,t) = tempF0(:,:,t)-mi;
% end

ma = max(tempF0,[],3);  %Do this first so that it has dF/F units (for later)
mi = min(tempF0,[],3);
%mapmag = (ma-mi)./(ma+mi);
tfmapmag = ma-mi;
tfmapmag = ma;

tfDomTens = zeros(size(tempF0));
for i = 1:length(tdom)
    tfDomTens(:,:,i) = log2(tdom(i));
end

[mi dum] = min(tempF0,[],3);
for i = 1:length(tdom)
    tempF0(:,:,i) = phi(tempF0(:,:,i)-mi);
end

sumMat = sum(tempF0,3);
for i = 1:length(tdom)
    tempF0(:,:,i) = tempF0(:,:,i)./sumMat;
end

tfMap = sum(tempF0.*tfDomTens,3); %center of mass
tfMap = 2.^tfMap;  %put back in the right units
tfMap = stim_framerate./tfMap; %cyc/sec.

%tfMap = (tempF0(:,:,1)./tempF0(:,:,4));

%tfMap = (tempF0(:,:,1)-tempF0(:,:,4))./(tempF0(:,:,1)+tempF0(:,:,4));
%tfMap2 = (tempF0(:,:,2)-tempF0(:,:,3))./(tempF0(:,:,2)+tempF0(:,:,3));
% a = tempF0(:,:,1)+ tempF0(:,:,2);
% b = tempF0(:,:,3)+ tempF0(:,:,4);
%tfMap = (a-b)./(a+b);

% tfMap(find(tfMap<-1)) = -1;
% tfMap(find(tfMap>1)) = 1;
    

%%
areaNames = varargin{1};
for i = 1:length(areaNames)
    if strcmp(areaNames{i},'V1')
        V1id = i;
    end
end

imlabel = bwlabel(imstate.areaBounds,4);
areaID = unique(imlabel);
V1bw = zeros(size(imstate.areaBounds));
V1bw(find(imlabel == V1id)) = 1;

x = 1:size(V1bw,2);
y = 1:size(V1bw,1);

CoMX = round(sum(x.*sum(V1bw)/sum(V1bw(:))));
CoMY = round(sum(y'.*sum(V1bw,2)/sum(V1bw(:))));
% 
imstate.fmaps{1} = imstate.fmaps{1} - imstate.fmaps{1}(CoMY,CoMX);
imstate.fmaps{2} = imstate.fmaps{2} - imstate.fmaps{2}(CoMY,CoMX);


%% Set Bounds

%First make a mask that only includes the areas in 'areaOrder'
areaOrder{1} = 'V1';
areaOrder{2} = 'PM';
areaOrder{3} = 'LM';
areaOrder{4} = 'AL';
areaOrder{5} = 'RL';

imlabel = bwlabel(imstate.areaBounds,4);
aBound = zeros(size(imstate.fmaps{1}));
for j = 1:length(areaOrder)
    
    aID = [];
    for i = 1:length(areaNames)
        if strcmp(areaNames{i},areaOrder{j})
            aID = i;
        end
    end    
    if ~isempty(aID)
        id = find(imlabel == aID);
    end
    aBound(id) = 1;
    
end

%Now remove pixels based on certain criteria

%Limit eccentricity
EmaxH = 80;
EmaxV = 50;
id = find(imstate.fmaps{1}>EmaxH | imstate.fmaps{1}<-EmaxH | imstate.fmaps{2}>EmaxV | imstate.fmaps{2}<-EmaxV);
aBound(id) = 0;

Emax = 50;
%Remove pixels with pad parameter fits
id = find(isnan(colorMap));
aBound(id) = 0;

%Smooth the masks edges
% SE = strel('disk',4,0);
% aBound = imopen(aBound,SE);
% SE = strel('disk',2,0);
% aBound = imclose(aBound,SE);

close

plotImages(colorMap,tfMap,cdom,sort(stim_framerate./tdom),imstate,aBound,Emax)

areaScatter(colorMap,colorF0,tfMap,cdom,tdom,imstate,aBound,areaOrder,areaNames,Emax)

%% Plot images

function plotImages(colorMap,tfMap,cdom,tdom,imstate,aBound,Emax)

global pixpermm

xdom = (0:(size(imstate.fmaps{1},2)-1))/pixpermm;
ydom = (0:(size(imstate.fmaps{1},1)-1))/pixpermm;

figure, 
ax1 = subplot(1,4,1);
imagesc(xdom,ydom,colorMap,'AlphaData',aBound), title('Direction in S/M plane'), axis image, colorbar
%plotExLocs((E_yx-1)/pixpermm)
hold on 
contour(xdom,ydom,imstate.areaBounds,[.5 .5]*180,'k')
caxis([0 180])
colormap(ax1,hsv)

plotRange = log2([tdom(1) tdom(end)]); %octaves
ax2 = subplot(1,4,2);
imagesc(xdom,ydom,log2(tfMap),'AlphaData',aBound,plotRange), title('temporal frequency'), axis image, colorbar
% plotExLocs((E_yx-1)/pixpermm)
hold on
contour(xdom,ydom,imstate.areaBounds,[.5 .5],'k')

clear domcell


iddom = linspace(plotRange(1),plotRange(2),5);
for i = 1:length(iddom)
    domcell{i} = num2str(round(iddom(i)*10)/10); %the label
end
colorbar(ax2,'YTick',iddom,'YTickLabel',domcell)

% caxis([min(pS(:)) max(pS(:))])
% Smap = [64 17 64]/64; Smap = (0:63)'*Smap/63;
colormap(ax2,jet)

%bnd = imstate.areaBounds*(max(pS(:))-min(pS(:))) + min(pS(:));
%contour(xdom,ydom,bnd,[.5 .5]*(max(pS(:))-min(pS(:))),'k')

ax3 = subplot(1,4,3);
imagesc(xdom,ydom,imstate.fmaps{1},'AlphaData',aBound), title('hor ret'), axis image, axis image, colorbar
%plotExLocs((E_yx-1)/pixpermm)
hold on 
contour(xdom,ydom,imstate.areaBounds*10,[5 5],'k')
caxis([-Emax Emax])
colormap(ax3,jet)

ax4 = subplot(1,4,4);
imagesc(xdom,ydom,imstate.fmaps{2},'AlphaData',aBound, [-Emax Emax]), title('vert ret'), axis image, axis image, colorbar
%plotExLocs((E_yx-1)/pixpermm)
hold on 
contour(xdom,ydom,imstate.areaBounds*10,[5 5],'k')
%hold on 
%contour(xdom,ydom,imstate.fmaps{1}.*aBound,[-40:20:40],'k','LineWidth',.01)
caxis([-Emax Emax])
colormap(ax4,jet)


%% Scatter plots of color vs. vertical retinotopy, for each area

function areaScatter(colorMap,colorF0,tfMap,cdom,tdom,imstate,aBound,areaOrder,areaNames,Emax)

%vertBinEdges = [-45:15:45];

mapplot = abs(sin(colorMap*pi/180));

mapplot = colorF0(:,:,3)./(colorF0(:,:,1)+colorF0(:,:,3));

imlabel = bwlabel(imstate.areaBounds,4);
areaID = unique(imlabel);

areaOrder{1} = 'V1';
areaOrder{2} = 'PM';
areaOrder{3} = 'LM';
areaOrder{4} = 'AL';
areaOrder{5} = 'RL';


clear sz
for j = 1:length(areaOrder)
    
    aID = [];
    for i = 1:length(areaNames)
        if strcmp(areaNames{i},areaOrder{j})
            aID = i;
        end
    end
    
    if ~isempty(aID)
        
        id = find(imlabel == aID & abs(imstate.fmaps{1})<Emax & abs(imstate.fmaps{2})<Emax);
        vertdum = imstate.fmaps{2}(id);
        cmapdum = mapplot(id);
        %         id = find(cmapdum>-.1 & cmapdum<1.5);
        %         vertdum = vertdum(id);
        %         cmapdum = cmapdum(id);
        
        id = find(isnan(cmapdum));
        cmapdum(id) = [];
        vertdum(id) = [];
        
        figure(69)
        if j == 1
            subplot(2,3,[1 4])
        elseif j<4
            subplot(2,3,j)
        elseif j>3
            subplot(2,3,j+1)
        end
        %subplot(5,1,j)
        
        DataS.vert_vec{aID} = vertdum;
        DataS.cmap_vec{aID} = cmapdum;
        DataS.cmap = colorMap;
        DataS.vertmap = imstate.fmaps{2};
        DataS.hormap = imstate.fmaps{1};
        DataS.areaBounds = imstate.areaBounds;
        
        %Plot density map
        [mat xdom ydom] = smoothscatter(vertdum,cmapdum,1,.05,[-Emax Emax],[-.5 1.5]);
        
        mat = mat-min(mat(:));
        mat = mat/max(mat(:));
        imagesc(xdom,ydom,1-mat), colormap gray
        xlim([-Emax Emax]),
        ylim([-.1 1.1])
        axis xy
        
        %Plot sigmoidal fit
        if j == 1 %do it only for V1
            [param cmapHat varaccount] = Sigfit(vertdum,cmapdum);
            sigdom = -Emax:Emax;
            d = sigdom-param(1);
            sigffit = param(3)*(1./(1 + exp(-d*param(2))) - .5) + param(4);
            DataS.sigfit{aID} = param;
        end
        hold on, plot(sigdom,sigffit,'b')
        
        slopeAtorigin = round(param(3)*param(2)/4*1000)/10; %perc/deg
        
        %Plot linear fit
        H = [vertdum ones(size(vertdum))];
        y = cmapdum;
        xhat = inv(H'*H)*H'*y;
        dom = -Emax:Emax;
        lfit = dom*xhat(1) + xhat(2);
        
        H = [cmapdum ones(size(cmapdum))];
        y = vertdum;
        xhatINV = inv(H'*H)*H'*y;
        xhat2(1) = 1/xhatINV(1);
        xhat2(2) = -xhatINV(2)/xhatINV(1);
        lfit2 = dom*xhat2(1) + xhat2(2);
        
        xhat = (xhat+xhat2')/2;
        
        linfit = dom*xhat(1) + xhat(2);
        DataS.linfit{aID} = xhat;
        
        hold on, plot(dom,linfit,'y')
        
        SperDeg = round(xhat(1)*1000)/10; %percent per degree
        
        [r p] = corrcoef(vertdum,cmapdum);
        r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
        
        titlestr = [areaOrder{j} '; r=' num2str(r) '; p=' num2str(p) '; %/deg=' num2str(SperDeg)]
        if j == 1
            titlestr = [titlestr '; ' num2str(slopeAtorigin)] ;
        end
        title(titlestr);
        %title([areaOrder{j} '; %var=' num2str(round(varaccount*100)/100)]);
        
        RetinaS = getRetinaGradient(-Emax:Emax);
        hold on,
        plot(-Emax:Emax,RetinaS/100)
        
        
        %Create bins based on percentiles
        vertBinEdges = [ ];
        prcDom = linspace(0,100,8)
        for i = 1:length(prcDom)
            vertBinEdges = [vertBinEdges prctile(vertdum,prcDom(i))];
        end
        %vertBinEdges = linspace(-Emax,Emax,8);
        
        clear mapmu mapsig vertDomain
        for k = 1:length(vertBinEdges)-1
            
            id = find(vertdum>=vertBinEdges(k) & vertdum<vertBinEdges(k+1));
            cmapdumdum = cmapdum(id);
            vertdumdum = vertdum(id);
            
            mi = prctile(cmapdumdum,5);
            ma = prctile(cmapdumdum,95);
            idOutlier = find(cmapdumdum<mi | cmapdumdum>ma);
            vertdumdum(idOutlier) = [];
            cmapdumdum(idOutlier) = [];
            
            vertDomain(k) = nanmedian(vertdumdum);
            mapmu(k) = nanmedian(cmapdumdum);
            mapsig(k) = nanstd(cmapdumdum);
            
        end
        
        hold on,
        errorbar(vertDomain,mapmu,mapsig,'r')
        
        ylabel('%S')
        xlabel('vertical eccentricity')
        
        ma = max(mapmu)+2*max(mapsig);
        mi = min(mapmu)-2*max(mapsig);
        
        axis square
        
        
        
        %ylim([mi ma])
        %ylim([-3 3])
        %xlim([-Emax Emax])
        %ylim([0 90])
        
        figure(70)
        subplot(3,2,j+1)
        
        plot([0 45],[-.0 -.0],'k')
        hold on
        plot([-45 -45],[0 1],'k')
        hold on
        plot([0 0],[0 0],'ok','MarkerSize',5)
        hold on
        
        hold on, plot(sigdom,sigffit,'b')
        
        errorbar(vertDomain,mapmu,mapsig,'r')
        xlim([-Emax Emax]),
        ylim([-.1 1.1])
        axis square
        title(titlestr);
        axis off
        
    end
end

%% Plot examples

function plotExamples(imstate,imMatF,E_yx)

HVpos = [0 -20; 0 0; 0 20];
HVpos = [0 -15; 0 15]

for i = 1:length(varargin{1})
    if strcmp(varargin{1}{i},'V1')
        V1id = i;
        break
    end
end

imlabel = bwlabel(imstate.areaBounds,4);
id = find(imlabel == V1id & aBound == 1);
V1bw = NaN*zeros(size(imstate.areaBounds));
V1bw(id) = 1;

sig = 10; %Smooth before identifing location of example
hh = fspecial('gaussian',size(imMat{1,1}),sig);
Hdum = ifft2(abs(fft2(hh)).*fft2(imstate.fmaps{1}));
Vdum = ifft2(abs(fft2(hh)).*fft2(imstate.fmaps{2}));
Hdum = V1bw.*Hdum;
Vdum = V1bw.*Vdum;
E_yx = [];


if isempty(E_yx)

    for i = 1:size(HVpos,1)
        D = sqrt((Vdum-HVpos(i,2)).^2 + (Hdum-HVpos(i,1)).^2);
        mi = min(D(:));
        [E_yx(i,1) E_yx(i,2)] = find(D == mi);
    end
    
end

%cid = 2:7;
%gdomH = cdom(cid);
%y = zeros(length(tdom)*length(gdomH),size(imMatF{1},1), size(imMatF{1},2));

k = 1;
clear y
for c = 1:length(cdom)
    for t = 1:length(tdom)
        y(k,:,:) = imMatF{t}(:,:,c); %layered cake
        k = k+1;
    end
end
Npts = k-1;

if normflag
    
   ynorm = sqrt(sum(y.^2));
   for i = 1:size(y,1)
       y(i,:,:) = y(i,:,:)./ynorm;
   end
    
end

figure
for i = 1:size(E_yx,1)
    E_vecloc = E_yx(i,1)+(E_yx(i,2)-1)*size(imstate.areaBounds,1);
    kern = squeeze(y(:,E_yx(i,1),E_yx(i,2)));
    kern = reshape(kern,length(tdom),length(cdom));
    kernfit = yhat(:,E_vecloc);
    kernfit = reshape(kernfit,length(tdom),length(cdom));
    ma = max([kern(:); kernfit(:)]);
    mi = min([kern(:); kernfit(:)]);
    
    subplot(4,size(E_yx,1),i)
    imagesc(kern,[mi ma]), axis xy
    title(['H=' num2str(HVpos(i,1)) '; V=' num2str(HVpos(i,2))])
    subplot(4,size(E_yx,1),size(E_yx,1)+i)
    imagesc(kernfit,[mi ma]), axis xy
    title(['M=' num2str(pM(E_vecloc)) '; S=' num2str(pS(E_vecloc))]) 
    colormap parula
    
    colorarray = parula;
    subplot(4,size(E_yx,1),2*size(E_yx,1)+i)
    plot([0 0],[-1 1],'--k'), hold on, plot([-1 1],[0 0],'--k'), hold on
    for j = 1:length(gdomM(:))
        kernval = kern(j);
        analIDX = (kernval-mi)/(ma-mi);
        discIDX = round(analIDX*(length(colorarray(:,1))-1))+1;
        
        plot((gdomM(j)),(bdomM(j)),'.','MarkerSize',20,'MarkerEdgeColor',colorarray(discIDX,:)), xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
        hold on        
        
        plot((-gdomM(j)),(-bdomM(j)),'.','MarkerSize',20,'MarkerEdgeColor',colorarray(discIDX,:)), xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
        hold on  
        
    end
    axis square
    xlim([-1.1 1.1])
    ylim([-1.1 1.1])
    axis off
    
    subplot(4,size(E_yx,1),3*size(E_yx,1)+i)
    [Mdomdum Sdomdum] = meshgrid(-1:.01:1,-1:.01:1);
    fitIm = abs(Mdomdum)*pM(E_vecloc) + abs(Sdomdum)*pS(E_vecloc);
    %fitIm = ((Mdomdum.^2)*pM(E_vecloc) + (Sdomdum.^2)*pS(E_vecloc) + (Sdomdum.*Mdomdum)*colorMap(E_vecloc));
    %fitIm = sqrt(fitIm);
    imagesc(Mdomdum(1,:),Sdomdum(:,1)',fitIm)
    xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
    axis square
    hold on, plot([0 0],[-1 1],'--k'), hold on, plot([-1 1],[0 0],'--k'), hold on
    rat = pS(E_vecloc)/(pS(E_vecloc)+pM(E_vecloc));
    title(['%S = ' num2str(round(100*rat)/100)])
    axis off

end
colorbar

close 


%% Response matrix averaged for different eccentricities
% imlabel = bwlabel(imstate.areaBounds,4);
% areaID = unique(imlabel);
% 
% vertBinEdges = [-67.5:45:67.5];
% 
% clear gmean gstd
% for i = 2:length(areaID) %loop each area
%     for j = 1:length(vertBinEdges)-1 %Loop each bin of vertical retinotopy
%         id = find(imlabel == areaID(i) & imstate.fmaps{2}>vertBinEdges(j) & imstate.fmaps{2}<vertBinEdges(j+1));
%         for k = 1:length(tdom) %Loop each t_period
%             for l = 1:length(cdom) %Loop each theta
%                 dum = imMatF{k}(:,:,l);
%                 gmean{i-1}(j,k,l) = mean(dum(id));
%                 gstd{i-1}(j,k,l) = std(dum(id));
%             end
%         end
%     end
% end
% 
% for i = 1:length(vertBinEdges)-1    
%    leg{i} = [num2str(vertBinEdges(i)) 'to'   num2str(vertBinEdges(i+1))];
% end
% 
% retdom = (vertBinEdges(1:end-1)+vertBinEdges(2:end))/2;
% 
% figure,
% q = 1;
% for k = 1:length(tdom) %Loop each t_period
%     for i = 2:length(areaID) %loop each area
%         
%         subplot(length(tdom),length(areaID)-1,q)
%         mudum = squeeze(gmean{i-1}(:,k,:))';      
%         stddum = squeeze(gstd{i-1}(:,k,:))'; 
%         %imagesc(cdom,retdom,mudum'), xlabel('theta'), ylabel('retinotopy')
%         errorbar(cdom'*ones(1,size(mudum,2))./tdom(k),mudum,stddum), xlabel('col/tp')
%         
%         areastr = '';
%         if ~isempty(varargin)
%             areaNames = varargin{1};
%             areastr = areaNames{i-1};
%         end
%         
%          title([areastr '; tp = ' num2str(tdom(k))])
%          if q == 1
%              legend(leg)
%          end
%          %xlim([cdom(1)-.1 cdom(end)+.1])
%         xlim([-5 5])
%         q = q+1;
%     end
% end
%%
% figure,
% q = 1;
% [gdomM bdomM] = meshgrid(cdom,tdom);
% ratDom = gdomM./bdomM;
% [ratDom,ratDomID] = sort(ratDom(:));
% for i = 2:length(areaID) %loop each area
%     
%     subplot(1,length(areaID)-1,q)
%     mudum = squeeze(gmean{i-1}(:,:,:));  %ret x tp x col
%     stddum = squeeze(gstd{i-1}(:,:,:));
%     for j = 1:length(retdom)
%         mudumRet = squeeze(mudum(j,:,:)); %tp x col
%         stddumRet = squeeze(stddum(j,:,:));
%         %imagesc(cdom,retdom,mudum'), xlabel('theta'), ylabel('retinotopy')
%         %errorbar(cdom'*ones(1,size(mudum,2)),mudum,stddum), xlabel('theta')
%         mudumRet = mudumRet;
%         plot(ratDom,mudumRet(ratDomID))
%         hold on
%         %title(['tp = ' num2str(tdom(k)) '; area' num2str(i-1)])
%         %xlim([cdom(1)-.1 cdom(end)+.1])
%     end
%     q = q+1;
% end

%%
% figure,
% q = 1;
% for i = 2:length(areaID) %loop each area
%     for j = 1:length(vertBinEdges)-1 %Loop each eccentricity
%     
%         
%         subplot(length(vertBinEdges)-1,length(areaID)-1,q)
%         dum = squeeze(gmean{i-1}(j,:,:));        
%         imagesc(cdom,tdom,dum), ylabel('t_period'), xlabel('theta')
%          %title(['tp = ' num2str(tdom(k)) '; area' num2str(i-1)])
%          axis xy
%          
%          
%          %plot(cdom,dum')
%         
%         q = q+1;
%     end
% end
% 
% 
% %%
% hall = [];
% %for i = 1:length(imMatF) %Loop each blue gain
% for i = 1:3 %Loop each blue gain
%     figure
%     [dum id] = min(imMatF{i},[],3);
%     miLoc = cdom(id)/tdom(i);
%     for j = 1:length(vertBinEdges)-1 %Loop each eccentricity
%         
%         id = find(imlabel == areaID(i) & imstate.fmaps{2}>vertBinEdges(j) & imstate.fmaps{2}<vertBinEdges(j+1));
%         h = miLoc(id);
%         hall = [hall; h(:)];
%         %subplot(1,length(vertBinEdges)-1,j),hist(h,7)
%         
%     end
% end


function val = getVal(pstring,cond)

global Analyzer

Nparam = length(Analyzer.loops.conds{1}.val);

for k = 1:Nparam
    if strcmp(Analyzer.loops.conds{cond}.symbol{k},pstring)
        val = Analyzer.loops.conds{cond}.val{k};
        break
    end
end

function plotExLocs(E_yx)

for i = 1:size(E_yx,1)
    hold on
    plot(E_yx(i,2),E_yx(i,1),'ok')
end