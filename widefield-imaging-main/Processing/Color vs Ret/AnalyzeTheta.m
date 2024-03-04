function DataS = AnalyzeTheta(imstate,f0m,E_yx,Gtrx,varargin)

global Analyzer pixpermm

odom = getdomain('ori');
cdom = getdomain('theta');
%tdom = getdomain('t_period');
tdom = 1;

imstate.areaBounds(find(imstate.areaBounds<1)) = 0;
imstate.areaBounds(find(~imstate.bw)) = 0;

%Smooth the boundary edges
sig = 3;
hh = fspecial('gaussian',size(imstate.areaBounds),sig);
imstate.areaBounds = ifft2(fft2(imstate.areaBounds).*abs(fft2(hh)));
imstate.areaBounds = (sign(imstate.areaBounds-.5)+1)/2;

DataS.areaNames = varargin{1};

pixpermm = 125;

%% Build cell matrix of responses to each color

imMat = cell(length(tdom),length(cdom));
for cond = 1:length(f0m)-1  
    
    theta = getVal('theta',cond);
    cid = find(cdom == theta);
    
    %t_period = getVal('t_period',cond);
    %tid = find(tdom == t_period);
    tid = 1;
    
    if isempty(imMat{tid,cid})
        imMat{tid,cid} = 0;
    end
    
    imMat{tid,cid} = imMat{tid,cid}+f0m{cond}/length(odom);
    
    col(tid,cid) = theta;
    %tp(tid,cid) = t_period;
    
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

colorMapdum = colorF0(:,:,2)./(colorF0(:,:,1)+colorF0(:,:,2));
%colorMapdum = real(log2(colorF0(:,:,2)./(colorF0(:,:,1))));

clims = [-0 1]

if ~isempty(Gtrx)
    Ref_Im = imref2d(size(imstate.fmaps{1}));
    colorMap = imwarp(colorMapdum,Gtrx,'OutputView',Ref_Im);
    id = find(colorMap==0);  %Values at 0 are the edges of the image after the transformation.
    colorMap(id) = NaN;
else
    colorMap = colorMapdum;
end

id = find(colorMap< clims(1)-.1);
colorMap(id) = NaN;
id = find(colorMap> clims(2)+.1);
colorMap(id) = NaN;

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
EmaxH = 60;
EmaxV = 60;
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

plotImages(colorMap,imstate,aBound,Emax,clims)

DataS = areaScatter(colorMap,imstate,aBound,areaOrder,areaNames,Emax)

%% Plot images

function plotImages(colorMap,imstate,aBound,Emax,clims)

global pixpermm

xdom = (0:(size(imstate.fmaps{1},2)-1))/pixpermm;
ydom = (0:(size(imstate.fmaps{1},1)-1))/pixpermm;

figure, 
ax1 = subplot(1,3,1);
imagesc(xdom,ydom,colorMap,'AlphaData',aBound,[clims(1) clims(2)]), title('%S opsin'), axis image, colorbar
%plotExLocs((E_yx-1)/pixpermm)
hold on 
contour(xdom,ydom,imstate.areaBounds,[.5 .5]*180,'k')
caxis([0 1])
colormap(ax1,jet)

clear domcell

% caxis([min(pS(:)) max(pS(:))])
% Smap = [64 17 64]/64; Smap = (0:63)'*Smap/63;

%bnd = imstate.areaBounds*(max(pS(:))-min(pS(:))) + min(pS(:));
%contour(xdom,ydom,bnd,[.5 .5]*(max(pS(:))-min(pS(:))),'k')

ax3 = subplot(1,3,2);
imagesc(xdom,ydom,imstate.fmaps{1},'AlphaData',aBound), title('hor ret'), axis image, axis image, colorbar
%plotExLocs((E_yx-1)/pixpermm)
hold on 
contour(xdom,ydom,imstate.areaBounds*10,[5 5],'k')
caxis([-Emax Emax])
colormap(ax3,jet)

ax4 = subplot(1,3,3);
imagesc(xdom,ydom,imstate.fmaps{2},'AlphaData',aBound, [-Emax Emax]), title('vert ret'), axis image, axis image, colorbar
%plotExLocs((E_yx-1)/pixpermm)
hold on 
contour(xdom,ydom,imstate.areaBounds*10,[5 5],'k')
%hold on 
%contour(xdom,ydom,imstate.fmaps{1}.*aBound,[-40:20:40],'k','LineWidth',.01)
caxis([-Emax Emax])
colormap(ax4,jet)


%% Scatter plots of color vs. vertical retinotopy, for each area

function DataS = areaScatter(colorMap,imstate,aBound,areaOrder,areaNames,Emax)

%vertBinEdges = [-45:15:45];

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
        cmapdum = colorMap(id);
        %         id = find(cmapdum>-.1 & cmapdum<1.5);
        %         vertdum = vertdum(id);
        %         cmapdum = cmapdum(id);
        
        id = find(isnan(cmapdum));
        cmapdum(id) = [];
        vertdum(id) = [];
        
        figure
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
            [param cmapHat varaccount] = Sigfit(vertdum,cmapdum,[0 .1 1 0]);
            sigdom = -Emax:Emax;
            d = sigdom-param(1);
            sigffit = param(3)./(1 + exp(-d*param(2))) + param(4);
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
            
            mi = prctile(cmapdumdum,2);
            ma = prctile(cmapdumdum,98);
            idOutlier = find(cmapdumdum<mi | cmapdumdum>ma);
            vertdumdum(idOutlier) = [];
            cmapdumdum(idOutlier) = [];
            
            vertDomain(k) = nanmedian(vertdumdum);
            mapmu(k) = nanmedian(cmapdumdum);
            mapsig(k) = nanstd(cmapdumdum);
            
        end
        
        DataS.ebarmu{aID} = mapmu;
        DataS.ebarsig{aID} = mapsig;
        DataS.ebardom{aID} = vertDomain;
        
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
        
        figure
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

