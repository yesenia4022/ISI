pF0_ISI

%%
global AUE ACQinfo Analyzer datadir

anim = 'xt7'; 
expt = 'u009_004'; %epi ret expt
ISIflag = 1;


anim = 'xt8'
expt = 'u003_009'; 
ISIflag = 1;
pixpermm = 79; 

anim = 'yd3'
expt = 'u001_002'; 
pixpermm = 79; 
ISIflag = 0;


%expt = 'u000_002'; %ori ocdom experiment (right hemi)

anadir = 'c:\AnalyzerFiles\'; %partial path for analyzer file
datadir = 'h:\ISIdata\'; %partial path for data files 

anadir = [anadir anim '\' ];
datadir = [datadir anim '\' expt '\'];  %Make path for data files

load([anadir anim '_' expt '.analyzer'],'-mat') %load analyzer file

%Get the syncInfo and put it in the Analyzer structure for convenience
i = 1;
while exist(['syncInfo' num2str(i)])    
   eval(['syncInfoAll{i} = syncInfo' num2str(i) ';']) 
   i = i+1;   
end
Analyzer.syncInfo = syncInfoAll;


AUE = [anim '_' expt]; %loaded unit experiment 'u000_000'
setacqinfo(1)  %Set the global ACQinfo structure (contains trial info

%% Set trial averaging window


%ISI
if ISIflag
    baselineWin = [-500 500]; %ms [start stop]
    ResponseWin = [1000 getparam('stim_time')*1000+2000];

%GCaMP
else
    baselineWin = [-500 0]; %ms [start stop]
    ResponseWin = [100 getparam('stim_time')*1000+500]; 
end


B_Flim = getframeidx(baselineWin,1); %Identify frame number delimits for baseline in each trial
T_Flim = getframeidx(ResponseWin,1); %Identify frame number delimits for response in each trial


%% Select region-of-interest

dum = GetTrialData([1000 1000],1); %Get a single frame
figure,imagesc(dum), colormap gray
bw = roipoly;


%% Get a baseline and response image for each trial
N = ACQinfo.linesPerFrame;
N(2) = ACQinfo.pixelsPerLine;
RMat = zeros(N(1),N(2),getnotrials);
BMat = zeros(N(1),N(2),getnotrials);
clear tc Rsep
for tno = 1:getnotrials
    tno
    Rxt = GetTrialData([-inf inf],tno);        %Get the entire trial
    
    if ISIflag %Invert, b/c ISI is a negative signal.
        bitdepth = 16;  %panda
        Rxt = 2^bitdepth - Rxt;
    end
    
    %Get time course for the trial
%     for t = 1:size(Rxt,3)
%         dum = Rxt(:,:,t);
%         tc(t,tno) = mean(dum(bw));
%     end

%     tfit = polyfitter(tc(:,tno)',1);
%     for t = 1:size(Rxt,3)
%         Rxt(:,:,t) = Rxt(:,:,t)-tfit(t) + mean(tfit);
%     end
  
        
    RMat(:,:,tno) = mean(Rxt(:,:,T_Flim(1):T_Flim(2)),3); %Store the mean response image of this trial
    BMat(:,:,tno) = mean(Rxt(:,:,B_Flim(1):B_Flim(2)),3); %Store the mean baseline image of this trial
    
    
end

%% ID saturated pixels and remove from ROI
% dum = max(RMat,[],3);
% idSat = find(dum>(2^bidDepth-20));
% bwSat = ones(N,N);
% bwSat(idSat) = 0;
% bw = bw.*bwSat;
% figure,imagesc(bw)

%% De-trend global signal
% clear Rsep
% for tno = 1:size(RMat,3)
%     dum = RMat(:,:,tno);
%     Rsep(:,tno) = dum(bw);
%     %Rsep(:,tno) = dum(:);
% end
% 
% %mu = prctile(Rsep(:),0);
% %[u s v] = svd(phi(Rsep-mu),'econ');
% 
% %Rx = u(:,1)*sign(v(1,1));
% %Rt = v(:,1)'*sign(u(1,1));
% Rt = sqrt(mean(Rsep.^2));
% Rtfit = polyfitter(Rt,5);
% RMat2 = RMat;
% BMat2 = BMat;
% for t = 1:size(RMat,3)
%     RMat2(:,:,t) = RMat2(:,:,t)/Rt(t);
%     BMat2(:,:,t) = BMat2(:,:,t)/Rt(t);
% end

%% Or, just subtract baseline on each trial?
% 
% RMat2 = RMat2;
%RMat2 = (RMat-BMat) + mean(BMat,3);
RMat2 = (RMat-BMat)./BMat;

% RBlank_Bsub = (Rblank3-Bblank)./BMat;

%RMat_Bsub = RMat_Bsub;
%RBlank_Bsub = (Rblank);



%% Organize by condition number

RMat3 = zeros(N(1),N(2),getnoconditions);
BMat3 = zeros(N(1),N(2),getnoconditions);

for c = 1:getnoconditions
    c
    nr = getnorepeats(c);
  
    clear tnos
    for r = 1:nr
        tnos(r) = Analyzer.loops.conds{c}.repeats{r}.trialno;
    end
    RMat3(:,:,c) = nanmean(RMat2(:,:,tnos),3); %mean over repeats
    BMat3(:,:,c) = nanmean(BMat(:,:,tnos),3);
    
end

blankflag = stimblank(getnoconditions);

% if blankflag
%     Rblank = RMat(:,:,end); %response to blanks
%     RMat(:,:,end) = [];
%     
%     Bblank = BMat(:,end); %baseline to blanks
%     BMat(:,:,end) = [];
% end


%%  All the code up until here could be applied to any experiment.  
% Here is where it starts to assume certain looper variable 

%get the index within the looper
for p = 1:length(Analyzer.L.param)
    if strcmp('s_phase2',Analyzer.L.param{p}{1})
        phaseDim = p;
    end
    if strcmp('ori2',Analyzer.L.param{p}{1})
        oriDim = p;
    end
    
end


%valmat is a matrix containing the looper values for each condition
clear valMat
blankflag = stimblank(getnoconditions);
for c = 1:(getnoconditions-blankflag)   
    for p = 1:length(Analyzer.loops.conds{c}.val)
        valMat(p,c) = Analyzer.loops.conds{c}.val{p};
    end    
end

oridom = unique(valMat(oriDim,:));
phasedom = unique(valMat(phaseDim,:));


%% Get ROI
%figure,imagesc(mean(RMat3,3))
%bw = roipoly;


%% Plot mean response to each condition


figure
k = 1;
for ori = 1:length(oridom)
    for p = 1:length(phasedom)
   
        id = find(valMat(oriDim,:) == oridom(ori)  & valMat(phaseDim,:) == phasedom(p));
        imdum = RMat3(:,:,id);
        
        subplot(length(oridom),length(phasedom),k)
        
        
        %imdum = imrotate(imdum.*bw,-30,'crop');      
        
        imagesc(imdum.*bw,[0 .1]), colorbar
        %plot(imdum(180,:))
        
        k = k+1;
        
    end
    
end

%% Make the vector map at each orientation

xdom_mm = (0:(N(2)-1))/pixpermm;
ydom_mm = (0:(N(1)-1))/pixpermm;

figure
k = 1;
for ori = 1:length(oridom)
   imVec{ori} = zeros(N(1),N(2));
   im0 = zeros(N(1),N(2));
    for p = 1:length(phasedom)
   
        id = find(valMat(oriDim,:) == oridom(ori)  & valMat(phaseDim,:) == phasedom(p));
        imdum = RMat3(:,:,id);
        
        imdum = imdum.*bw;
        
        imVec{ori} = imVec{ori} + imdum*exp(1i*phasedom(p)*pi/180);
        im0 = im0+imdum;
       
        k = k+1;
        
    end
    imVec{ori} = imVec{ori}./im0; %normalize
    %imVec(:,:,ori) = imVec(:,:,ori) - mean(imVecdum(find(bw)));
    
    imVec{ori}(bw == 0) = 0;
    
    %imVecdum = imVec{ori};
    %imVecdum = exp(1i*angle(imVecdum));
    %imVecdum = ROIfilt(imVecdum,bw,1,inf);
    %imVecdum = imVecdum-nanmean(imVecdum(find(bw)));
    
    
    mag = abs(imVec{ori});
    mi = prctile(mag(find(bw)),1);
    ma = prctile(mag(find(bw)),99);
    
    ax1 = subplot(length(oridom),2,2*ori-1);
    imagesc(xdom_mm, ydom_mm, abs(imVec{ori}),[0 .3]), colormap(ax1,'gray'), colorbar
    title('phase selectivity')
    
    ax2 = subplot(length(oridom),2,2*ori);
    imagesc(xdom_mm,ydom_mm,angle(imVec{ori})*180/pi / 360 / getparam('s_freq2'),[-180 180] / 360 / getparam('s_freq2')), colormap(ax2,'hsv'), colorbar
    title('relative bar position (deg of visual field)')
    axis image
    xlabel('mm')
end

%%  Get an ROI inside of V1 using plots above

%bw_TC = roipoly();

%You need to first click on the figure that you want to draw the ROI box.
%Right click "crop image".
figure
for ori = 1:length(oridom)
    subplot(2,1,ori)
    
    imVecdum = imVec{ori};
    
    imVecdum = exp(1i*angle(imVecdum));
    
    imVecdum = ROIfilt(imVecdum,bw,1,inf);
    
    
    
    imagesc(angle(imVecdum)), colormap hsv


end

title('pick ROI')

[dum rectROI] = imcrop();



%%
sig = 1;
%xmin ymin width height = rectROI
rectROI = round(rectROI);
x1 = rectROI(1);
x2 = rectROI(1)+rectROI(3)-1;
y1 = rectROI(2);
y2 = rectROI(2)+rectROI(4)-1;

Hret = exp(1i*angle(imVec{1}));
Hret = angle(ROIfilt(Hret,bw,sig,inf))*180/pi;

Vret = exp(1i*angle(imVec{2}));  %Turn them into unity vectors
Vret = angle(ROIfilt(Vret,bw,sig,inf))*180/pi;  %Smooth, then convert into angle again.

Hret = Hret(y1:y2,x1:x2);  %Crop
Vret = Vret(y1:y2,x1:x2);  

figure
subplot(1,2,1);
imagesc(Hret,[-180 180]), colormap hsv

subplot(1,2,2);
imagesc(Vret,[-180 180]), colormap hsv

%%

[xmesh ymesh] = meshgrid((x1:x2)-x1,(y1:y2)-y1);

[Hor_phase_mod Vert_phase_mod] = retinotopy_unwrap(Hret(:),Vret(:),[ymesh(:) xmesh(:)]);
Hor_phase_mod = reshape(Hor_phase_mod,size(Hret));
Vert_phase_mod = reshape(Vert_phase_mod,size(Vret));

Hor_mod = Hor_phase_mod/360/getparam('s_freq2'); %Convert to deg of visual field
Vert_mod = Vert_phase_mod/360/getparam('s_freq2');

%%
% [param ffit varacc] = Planefit(xmesh(:),ymesh(:),Hret(:)+180)
% ffit = reshape(ffit,size(Hret));
% figure,
% subplot(2,2,1)
% imagesc(ffit), colorbar
% subplot(2,2,2)
% imagesc(Hret+180), colorbar
% 
% [param ffit varacc] = Planefit(xmesh(:),ymesh(:),Vret(:)+180)
% ffit = reshape(ffit,size(Hret));
% subplot(2,2,3)
% imagesc(ffit), colorbar
% subplot(2,2,4)
% imagesc(Vret+180), colorbar, colormap hsv
% 
% 
% [param ffit] = planefitter_bruteforce(xmesh(:)/pixpermm,ymesh(:)/pixpermm,Vret(:)+180)

%%

figure

H = [xmesh(:)/pixpermm ymesh(:)/pixpermm ones(prod(size(xmesh)),1)]; %domain in mm
Hor_plane = inv(H'*H)*H'*Hor_mod(:); %parameters of the plane
Hor_hat = H*Hor_plane;
Vert_plane = inv(H'*H)*H'*Vert_mod(:);
Vert_hat = H*Vert_plane;

Hor_hat = reshape(Hor_hat,size(Hret));
Vert_hat = reshape(Vert_hat,size(Vret));

miH = min(Hor_hat(:));
maH = max(Hor_hat(:));
miV = min(Vert_hat(:));
maV = max(Vert_hat(:));

subplot(2,2,1);
imagesc(xmesh(1,:)/pixpermm,ymesh(:,1)/pixpermm,Hor_mod,[miH maH]), colormap jet
axis image
colorbar 

subplot(2,2,2);
imagesc(xmesh(1,:)/pixpermm,ymesh(:,1)/pixpermm,Vert_mod,[miV maV]), colormap jet
axis image
colorbar 

subplot(2,2,3);
imagesc(xmesh(1,:)/pixpermm,ymesh(:,1)/pixpermm,Hor_hat,[miH maH]), colormap jet
axis image
colorbar 

subplot(2,2,4);
imagesc(xmesh(1,:)/pixpermm,ymesh(:,1)/pixpermm,Vert_hat,[miV maV]), colormap jet
axis image
colorbar

%% Calculate magnification

dHdx = Hor_plane(1);
dHdy = Hor_plane(2);
dVdx = Vert_plane(1);
dVdy = Vert_plane(2);

Mag = sqrt(abs(dHdx*dVdy - dVdx*dHdy))  %deg/mm

Mag_hor = sqrt(dHdx^2 + dHdy^2) %deg/mm
Mag_vert = sqrt(dVdx^2 + dVdy^2)


%%

sig = .5;
clear imV1
for ori = 1:length(oridom)
   im0 = zeros(N(1),N(2));
    for p = 1:length(phasedom)
   
        id = find(valMat(oriDim,:) == oridom(ori)  & valMat(phaseDim,:) == phasedom(p));
        imdum = RMat3(:,:,id);
        
        imdum = ROIfilt(imdum,bw,sig,inf);
        
        imdum = imdum(y1:y2,x1:x2);
        
        imV1{ori}(:,:,p) = imdum;
       
    end
end

%% Plot the phase tuning curve from a grid of superpixels
ngrid = 5;

xgrid = round(linspace(1,size(imV1{1},2),ngrid+1));
ygrid = round(linspace(1,size(imV1{1},1),ngrid+1));
figure
k = 1;
ori = 2;
for i = 1:length(xgrid)-1
    for j = 1:length(ygrid)-1
        
        dumTens = imV1{ori}(ygrid(j):ygrid(j+1),xgrid(i):xgrid(i+1),:);
        phaseTC = squeeze(mean(mean(dumTens,1),2));
        
        Pref{ori}(i,j) = angle(sum(phaseTC.*exp(1i*phasedom(:)*pi/180)))*180/pi;
        
        subplot(ngrid,ngrid,k)
        plot(phasedom,phaseTC)
        ylim([0 max(phaseTC)])
        k = k+1;
        
    end
end

figure,imagesc(Pref{ori}), colormap hsv

%%
clear phasetc_vertAll phasetc_horAll

np = length(phasedom);
dp = phasedom(2) - phasedom(1);
dim = size(Vert_hat);
phasetcAll = zeros(np,prod(dim));

npI = 18;
phasedomI = linspace(0,360,npI+1);
phasedomI = phasedomI(1:end-1);
dpI = phasedomI(2) - phasedomI(1);

k = 1;
for i = 1:dim(1) %loop through every pixel
    for j = 1:dim(2)
   
        phasetc_vert = squeeze(imV1{2}(i,j,:));
        phasetc_hor = squeeze(imV1{1}(i,j,:));
        
       phasetc_vert = [phasetc_vert; phasetc_vert(1)];
       phasetc_hor = [phasetc_hor; phasetc_hor(1)];
        
        phasetc_vertI = interp1([phasedom 360],phasetc_vert',phasedomI);
        phasetc_horI = interp1([phasedom 360],phasetc_hor',phasedomI);
        
        shift_vert = round(npI/2 - Vret(i,j)/dpI);
        shift_hor = round(npI/2 - Hret(i,j)/dpI);
        
%         [dum idma] = min(phasetc);
%         shift = round(np/2) - idma;
        
        phasetc_vertI = circshift(phasetc_vertI,[0 shift_vert]);
        phasetc_horI = circshift(phasetc_horI,[0 shift_hor]);
        
        phasetc_vertAll(:,k) = phasetc_vertI;
        phasetc_horAll(:,k) = phasetc_horI;
        k = k+1;
        
    end
end
%% Plot the average RF of each pixel
dom = phasedomI/360/getparam('s_freq2');
figure

subplot(1,2,1)
muRF = mean(phasetc_horAll,2);
plot(dom,muRF/max(muRF),'o-')
ylim([0 1.1])
xlabel('degrees of visual field')
[param ffit varacc sigma] = Gaussfit_alias(phasedomI(:)',muRF(:)'/max(muRF),1)
[param ffit varacc sigma] = Gaussfit(phasedomI(:)',muRF(:)'/max(muRF),1)
sig_hor_deg = param(2)/360/getparam('s_freq2')
hold on 
plot(dom,ffit)
sig_hor_mm = sig_hor_deg/Mag_hor
title(['hor RF; sig deg=' num2str(round(sig_hor_deg)) '; sig mm=' num2str(round(sig_hor_mm*100)/100)])

subplot(1,2,2)
muRF = mean(phasetc_vertAll,2);
plot(dom,muRF/max(muRF),'o-')
ylim([0 1.1])
title('vertical receptive field')
xlabel('degrees of visual field')
[param ffit varacc sigma] = Gaussfit_alias(phasedomI(:)',muRF(:)'/max(muRF),1)
[param ffit varacc sigma] = Gaussfit(phasedomI(:)',muRF(:)'/max(muRF),1)
sig_vert_deg = param(2)/360/getparam('s_freq2')
hold on 
plot(dom,ffit)
sig_vert_mm = sig_vert_deg/Mag_vert
title(['vert RF; sig deg=' num2str(round(sig_vert_deg)) '; sig mm=' num2str(round(sig_vert_mm*100)/100)])

%%  Calculate point-image size










%%



