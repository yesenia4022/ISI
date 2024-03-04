function DataS = AnalyzeBGmatrix(imstate,E_yx,varargin)

global f0m Analyzer

odom = getdomain('ori');
gdom = getdomain('greengain');
bdom = getdomain('bluegain');

imstate.areaBounds(find(imstate.areaBounds<1)) = 0;
imstate.areaBounds(find(~imstate.bw)) = 0;

DataS.areaNames = varargin{1};

%% Build cell matrix of responses to each color combination

imMat = cell(length(bdom),length(gdom));
for cond = 1:length(f0m)-1  
    
    greengain = getVal('greengain',cond);
    gid = find(gdom == greengain);
    
    bluegain = getVal('bluegain',cond);
    bid = find(bdom == bluegain);
    
    if isempty(imMat{bid,gid})
        imMat{bid,gid} = 0;
    end
    
    imMat{bid,gid} = imMat{bid,gid}+f0m{cond}/length(odom);
    
    ggain(bid,gid) = greengain;
    bgain(bid,gid) = bluegain;
    
end

%spatial Filter
sig = 5;
hh = fspecial('gaussian',size(imMat{1,1}),sig);
imMatF = cell(1,size(imMat,1));
for i = 1:size(imMat,1)
    for j = 1:size(imMat,2)        
        imMatF{i}(:,:,j) = ifft2(abs(fft2(hh)).*fft2(imMat{i,j}));        
    end
end

%% Get inverted center-of-mass of greengain response curve for each pixel and bluegain 
% 
% figure
% for i = 1:length(bdom) %Loop each value of bluegain
%     
%     imdum = -imMatF{i};
%     ma = min(imdum,[],3);
%     for j = 1:size(imdum,3)      
%         imdum(:,:,j) = imdum(:,:,j)-ma;        
%     end
%     imdum = phi(imdum);
%     
%     su = sum(imdum,3);
%     for j = 1:size(imdum,3)      
%         imdum(:,:,j) = imdum(:,:,j)./su;        
%     end
%    
%     
%     mag = mean(imMatF{i},3);
%     mi = prctile(mag(:),1);    
%     mag = mag-mi;
%     mag(find(mag<0)) = 0;
%     
%     ma = prctile(mag(:),50);
%     mag = mag/ma;
%     mag(find(mag>1)) = 1;  
%     
%     
%     CM = 0;
%     for j = 1:length(gdom)        
%         CM = CM + gdom(j)*imdum(:,:,j);        
%     end
%     
%     mi = prctile(CM(:),3);
%     ma = prctile(CM(:),97);
%     
%     subplot(1,3,i),imagesc(CM,'AlphaData',mag,[-.6 .6]), axis image, colorbar
%     title(['Gnotch; Bgain = ' num2str(bdom(i))])
%     colormap jet
%     hold on, contour(imstate.areaBounds,[.5 .5],'k')
% 
% end

%% Plot maps of from particular points of the response matrix
% 
% Miso = (imMatF{1}(:,:,7)+imMatF{1}(:,:,1))/2; %[b g] = [.1 1]
% Siso = imMatF{2}(:,:,2); %[b g] = [.3 -1]
% figure,
% subplot(1,3,1)
% imagesc((Miso-Siso)./(Miso+ Siso),[-.8 .8]), title('M-S'), colorbar, axis image
% hold on, contour(imstate.areaBounds,[.5 .5],'k')
% 
% 
% IsoLum = imMatF{3}(:,:,3); %[b g] = [1 -.25]
% Lum = imMatF{2}(:,:,5);  %[b g] = [.3 .25]
% subplot(1,3,2)
% imagesc((IsoLum-Lum)./(IsoLum+Lum),[-.8 .8]),  title('(L-M)-(L+M)'), colorbar, axis image
% hold on, contour(imstate.areaBounds,[.5 .5],'k')
% 
% 
% blue = imMatF{3}(:,:,4); %[b g] = [0.1 1]
% green = (imMatF{1}(:,:,1) + imMatF{1}(:,:,7))/2;  %[g b] = [.1 1/-1]
% subplot(1,3,3)
% imagesc((green-blue)./(blue+green),[-.8 .8]),  title('green-blue'), colorbar, axis image
% colormap jet
% hold on, contour(imstate.areaBounds,[.5 .5],'k')
% 
% 
% greenTC = 0;
% for i = 1:length(imMatF)
%     blueTC(:,:,i) = mean(imMatF{i},3); %[b g] = [0.1 1]
%     
%     greenTC = greenTC+imMatF{i}/length(imMatF);
% end
% 
% blueH = bdom(:);

%% Convert stimulus space to M and S
[gdomM bdomM] = getConeContrast(gdom,bdom); %This converts the vectors into a matrix
%[gdomM bdomM] = meshgrid(gdom,bdom);
[idb idg] = find(gdomM>=-10);
idv = find(gdomM>=-10);

% figure,
% subplot(1,3,1), 
% imagesc(1:length(gdom),1:length(bdom),bdomM,[-1 1]), colormap gray, colorbar, xlabel('green LED contrast'), ylabel('nUV LED contrast'), title('S-opsin contrast')
% set(gca,'Ytick',1:length(bdom),'YTicklabels',bdom,'Xtick',1:length(gdom),'XTicklabels',gdom)
% axis xy square
% subplot(1,3,2), 
% imagesc(1:length(gdom),1:length(bdom),gdomM,[-1 1]), colormap gray, colorbar, xlabel('green LED contrast'), ylabel('nUV LED contrast'), title('S-opsin contrast')
% set(gca,'Ytick',1:length(bdom),'YTicklabels',bdom,'Xtick',1:length(gdom),'XTicklabels',gdom)
% axis xy square

%% Plot locations in color plane
% figure 
% plot([0 0],[-1 1],'--r'), hold on, plot([-1 1],[0 0],'--r'), hold on
% plot((gdomM(:)),(bdomM(:)),'.','MarkerSize',30,'MarkerEdgeColor',[0 0 0]), xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
% 
% xlim([-1.1 1.1])
% ylim([-1.1 1.1])
% axis square

%% Fit plane at each pixel


k =1;
%gid = 2:7;
%gdomH = gdom(gid);
%y = zeros(length(bdom)*length(gdomH),size(imMatF{1},1), size(imMatF{1},2));
clear y
for g = 1:length(gdom)
    for b = 1:length(bdom)
        if ~isempty(find(b == idb & g == idg));
            y(k,:,:) = imMatF{b}(:,:,g); %layered cake
            k = k+1;
        end
    end
end
Npts = k-1;

y = reshape(y,[Npts size(imMatF{1},1)*size(imMatF{1},2)]);


%[gdomM bdomM] = meshgrid(gdomH,bdom);

%H = [gdomM(:) bdomM(:) ones(length(gdomM(:)),1)];
H = [gdomM(idv) bdomM(idv)];

%H = [abs(gdomM(idv)) abs(bdomM(idv))];
%H = [gdomM(idv).^2 bdomM(idv).^2];
%y = y.^2;


xhat = inv(H'*H)*H'*y;

yhat = H*xhat;
varacc = (var(y)-var(yhat - y))./var(y);

alpha = reshape(xhat(1,:),[size(imMatF{1},1), size(imMatF{1},2)]); 
beta = reshape(xhat(2,:),[size(imMatF{1},1), size(imMatF{1},2)]);
varacc = reshape(varacc,[size(imMatF{1},1), size(imMatF{1},2)]);

alpha(find(~imstate.areaBounds)) = NaN;
beta(find(~imstate.areaBounds)) = NaN;

alpha(find(varacc<.2)) = NaN;
beta(find(varacc<.2)) = NaN;

%alpha = abs(alpha); beta = abs(beta);

% alpha(find(alpha<=0)) = eps;
% alpha(find(beta<=0)) = eps;
% 
% alpha = sqrt(alpha);
% beta = sqrt(beta);



%% Plot examples

%gid = 2:7;
%gdomH = gdom(gid);
%y = zeros(length(bdom)*length(gdomH),size(imMatF{1},1), size(imMatF{1},2));
k = 1;
clear y
for g = 1:length(gdom)
    for b = 1:length(bdom)
        if ~isempty(find(b == idb & g == idg));
            y(k,:,:) = imMatF{b}(:,:,g); %layered cake
            k = k+1;
        end
    end
end
Npts = k-1;

figure
for i = 1:size(E_yx,1)
    E_vecloc = E_yx(i,1)+(E_yx(i,2)-1)*size(imstate.areaBounds,1);
    kern = squeeze(y(:,E_yx(i,1),E_yx(i,2)));
    kern = reshape(kern,length(bdom),length(gdom));
    kernfit = yhat(:,E_vecloc);
    kernfit = reshape(kernfit,length(bdom),length(gdom));
    ma = max([kern(:); kernfit(:)]);
    mi = min([kern(:); kernfit(:)]);
    
    subplot(4,size(E_yx,1),i)
    imagesc(kern,[mi ma]), axis xy
    subplot(4,size(E_yx,1),size(E_yx,1)+i)
    imagesc(kernfit,[mi ma]), axis xy
    title(['alpha=' num2str(alpha(E_vecloc)) '; beta=' num2str(beta(E_vecloc))]) 
    colormap parula
    
    colorarray = parula;
    subplot(4,size(E_yx,1),2*size(E_yx,1)+i)
    plot([0 0],[-1 1],'--r'), hold on, plot([-1 1],[0 0],'--k'), hold on
    for j = 1:length(gdomM(:))
        kernval = kern(j);
        analIDX = (kernval-mi)/(ma-mi);
        discIDX = round(analIDX*(length(colorarray(:,1))-1))+1;
        
        plot((gdomM(j)),(bdomM(j)),'.','MarkerSize',20,'MarkerEdgeColor',colorarray(discIDX,:)), xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
        hold on        
    end
    axis square
    xlim([-1.1 1.1])
    ylim([-1.1 1.1])
    
    subplot(4,size(E_yx,1),3*size(E_yx,1)+i)
    [Mdomdum Sdomdum] = meshgrid(-1:.01:1,-1:.01:1);
    fitIm = abs(Mdomdum)*alpha(E_vecloc) + abs(Sdomdum)*beta(E_vecloc);
    imagesc(Mdomdum(1,:),Sdomdum(:,1)',fitIm)
    xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
    axis square
    
    title(['alpha=' num2str(round(10^4*alpha(E_vecloc))/10^4) '; beta=' num2str(round(10^4*beta(E_vecloc))/10^4) ' 10^-^4']) 
    
%     subplot(4,size(E_yx,1),3*size(E_yx,1)+i)
%     plot([0 0],[-1 1],'--r'), hold on, plot([-1 1],[0 0],'--k'), hold on
%     for j = 1:length(gdomM(:))
%         kernval = kernfit(j);
%         analIDX = (kernval-min(kernfit(:)))/range(kernfit(:));
%         discIDX = round(analIDX*(length(colorarray(:,1))-1))+1;
%         
%         plot((gdomM(j)),(bdomM(j)),'.','MarkerSize',20,'MarkerEdgeColor',colorarray(discIDX,:)), xlabel('M-opsin contrast'), ylabel('S-opsin contrast')
%         hold on        
%     end
%     axis square
%     xlim([-1.1 1.1])
%     ylim([-1.1 1.1])
end

%% Plot images

figure, 
subplot(3,2,1),
imagesc(imstate.fmaps{1},'AlphaData',imstate.areaBounds), title('hor ret'), axis image, axis image, colorbar
plotExLocs(E_yx)

subplot(3,2,2),
imagesc(imstate.fmaps{2},'AlphaData',imstate.areaBounds), title('vert ret'), axis image, axis image, colorbar
plotExLocs(E_yx)

subplot(3,2,3),
imagesc(alpha,'AlphaData',imstate.areaBounds), title('M weight'), axis image
plotExLocs(E_yx)
subplot(3,2,4),
imagesc(beta,'AlphaData',imstate.areaBounds), title('S weight'), axis image
plotExLocs(E_yx)


%colorDir = atan2(beta,alpha)*180/pi;
colorDir = (beta-alpha)./(beta+alpha);
%colorDir = real(log(beta./alpha));

%colorDir = beta./(abs(alpha)+abs(beta)); colorDir(find(colorDir>1 | colorDir<0)) = NaN;

% colorMag = sqrt(alpha.^2 + beta.^2);
% mi = prctile(colorMag(:),1);
% colorMag = colorMag-mi;
% colorMag(find(colorMag<0)) = 0;
% ma = prctile(colorMag(:),60);
% colorMag = colorMag/ma;
% colorMag(find(colorMag>1)) = 1;
subplot(3,2,5)
imagesc(colorDir, 'AlphaData', imstate.areaBounds,[-1 1]), axis image, colorbar
title('(S-M)/(S+M)')
plotExLocs(E_yx)
%hold on, contour(imstate.areaBounds,[.5 .5],'k')

%% Scatter plots of color vs. vertical retinotopy, for each area

%vertBinEdges = [-45:15:45];

mapplot = colorDir;
imlabel = bwlabel(imstate.areaBounds,4);
areaID = unique(imlabel);



figure,
clear sz
for i = 2:length(areaID)
   id = find(imlabel == areaID(i) & abs(imstate.fmaps{1})<45 & abs(imstate.fmaps{2})<45);
   vertdum = imstate.fmaps{2}(id);
   cmapdum = mapplot(id);
   
   subplot(length(areaID)-1,1,i-1)

   
   DataS.vert_vec{i-1} = vertdum;
   DataS.cmap_vec{i-1} = cmapdum;
   DataS.cmap = colorDir;
   DataS.vertmap = imstate.fmaps{2};
   DataS.hormap = imstate.fmaps{1};
   DataS.areaBounds = imstate.areaBounds;
   
   [mat xdom ydom] = smoothscatter(vertdum,cmapdum,1,.05,[-45 45],[-1.5 1.5]);
   
   mat = mat-min(mat(:));
   mat = mat/max(mat(:));
   imagesc(xdom,ydom,1-mat), colormap gray
   xlim([-45 45]), ylim([-1.1 1.1])
   axis xy
   
   [r p] = corrcoef(vertdum,cmapdum);
   r = round(r(1,2)*100)/100; p = round(p(1,2)*100)/100;
   
   areastr = '';
   if ~isempty(varargin)
       areaNames = varargin{1};
       areastr = areaNames{i-1};
   end
   title([areastr '; r=' num2str(r) '; p=' num2str(p)]);
   
   %Create bins based on percentiles
   vertBinEdges = [ ];
   prcDom = 0:20:100;
   for i = 1:length(prcDom)
       vertBinEdges = [vertBinEdges prctile(vertdum,prcDom(i))];
   end
   
   clear mapmu mapsig vertDomain
   for j = 1:length(vertBinEdges)-1
       
       id = find(vertdum>=vertBinEdges(j) & vertdum<vertBinEdges(j+1));
       cmapdumdum = cmapdum(id);
       vertdumdum = vertdum(id);
       
       mi = prctile(cmapdumdum,5);
       ma = prctile(cmapdumdum,95);
       idOutlier = find(cmapdumdum<mi | cmapdumdum>ma);
       vertdumdum(idOutlier) = [];
       cmapdumdum(idOutlier) = [];
       
       vertDomain(j) = nanmedian(vertdumdum);
       mapmu(j) = nanmedian(cmapdumdum);
       mapsig(j) = nanstd(cmapdumdum);
       
   end
   
   hold on,
   
   errorbar(vertDomain,mapmu,mapsig,'r')
   ylabel('(S-M)/(S+M)')
   xlabel('vertical eccentricity')   

   ma = max(mapmu)+2*max(mapsig);
   mi = min(mapmu)-2*max(mapsig);

   %ylim([mi ma])
   %ylim([-3 3])
   %xlim([-45 45])
   %ylim([0 90])
   sz(i) = length(id); 
end

[dum id] = sort(sz);
V1id = id(2);
V1pix = find(imlabel == areaID(V1id));

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
%         for k = 1:length(bdom) %Loop each bluegain
%             for l = 1:length(gdom) %Loop each greengain
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
% for k = 1:length(bdom) %Loop each bluegain
%     for i = 2:length(areaID) %loop each area
%         
%         subplot(length(bdom),length(areaID)-1,q)
%         mudum = squeeze(gmean{i-1}(:,k,:))';      
%         stddum = squeeze(gstd{i-1}(:,k,:))'; 
%         %imagesc(gdom,retdom,mudum'), xlabel('greengain'), ylabel('retinotopy')
%         errorbar(gdom'*ones(1,size(mudum,2))./bdom(k),mudum,stddum), xlabel('ggain/bgain')
%         
%         areastr = '';
%         if ~isempty(varargin)
%             areaNames = varargin{1};
%             areastr = areaNames{i-1};
%         end
%         
%          title([areastr '; bgain = ' num2str(bdom(k))])
%          if q == 1
%              legend(leg)
%          end
%          %xlim([gdom(1)-.1 gdom(end)+.1])
%         xlim([-5 5])
%         q = q+1;
%     end
% end
%%
% figure,
% q = 1;
% [gdomM bdomM] = meshgrid(gdom,bdom);
% ratDom = gdomM./bdomM;
% [ratDom,ratDomID] = sort(ratDom(:));
% for i = 2:length(areaID) %loop each area
%     
%     subplot(1,length(areaID)-1,q)
%     mudum = squeeze(gmean{i-1}(:,:,:));  %ret x bgain x ggain
%     stddum = squeeze(gstd{i-1}(:,:,:));
%     for j = 1:length(retdom)
%         mudumRet = squeeze(mudum(j,:,:)); %bgain x ggain
%         stddumRet = squeeze(stddum(j,:,:));
%         %imagesc(gdom,retdom,mudum'), xlabel('greengain'), ylabel('retinotopy')
%         %errorbar(gdom'*ones(1,size(mudum,2)),mudum,stddum), xlabel('greengain')
%         mudumRet = mudumRet;
%         plot(ratDom,mudumRet(ratDomID))
%         hold on
%         %title(['bgain = ' num2str(bdom(k)) '; area' num2str(i-1)])
%         %xlim([gdom(1)-.1 gdom(end)+.1])
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
%         imagesc(gdom,bdom,dum), ylabel('bluegain'), xlabel('greengain')
%          %title(['bgain = ' num2str(bdom(k)) '; area' num2str(i-1)])
%          axis xy
%          
%          
%          %plot(gdom,dum')
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
%     miLoc = gdom(id)/bdom(i);
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