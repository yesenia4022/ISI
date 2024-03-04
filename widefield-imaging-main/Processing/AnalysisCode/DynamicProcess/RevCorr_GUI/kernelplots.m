function kernelplots(kerns)

%adapted from Gkernelplots2.  Kernels are provided as input, instead of
%taken from the 'cellS'

global ACQinfo Analyzer cellS maskS G_RChandles

%%%%

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
[nID] = getNeuronMask;

%%%%

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
try
    load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')
catch
    load(['e:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')
end
%load(['F:\neurostuff\log_files\' expt],'frate')
if ~exist('frate')
    frate = 60;
end

domains = getSeqInfo;

%%%%%%%%%%%%%%%%%%%
oridom = domains{1}.oridom;
sfdom = domains{1}.sfdom;
phasedom = domains{1}.phasedom;
colordom = domains{1}.colordom;

for i = 1:length(colordom)
    for p = 1:length(nID)
        kernC{i}{p} = kerns{p};
        kernSigC{i}{p} = kerns{p};
    end    
end

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string') ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with dtau spacing, and end at an estimate of kernDel(2)
taudom = (0:length(kernC{1}{1}(1,1,:))-1)*acqPeriod;

delayWin = [0 300];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

Ncell = length(nID);
NT = getnotrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)

kernsmooth = getSmoother([.1 .5 1 .5 .1],[.5 1 .5],[.1 1 .1],taudom,oridom,sfdom);

%Smoother before taking ori curve
orismoother = getSmoother([.1 .5 1 .5 .1],1,[.2 1 .2],taudom,oridom,sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother([.2 .5 1 .5 .2],[.3 1 .3],1,taudom,oridom,sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother(1,[.7 1 .7],[.4 1 .4],taudom,oridom,sfdom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorid = {'r','g','b'};
%plot sf curves

figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        kernplot = kernC{c}{p};
        %kernplot = squeeze(mean(kernplot,3)); %average over phase
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother)));        
        
        [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        
        tcsf = squeeze(kernplot(bestoriid,:,tau));
        
%         [dum sfpref(p)] = max(tcsf);
%         sfpref(p) = sfdom(sfpref(p));
%         sfmag(p) = (max(tcsf)-min(tcsf))/(max(tcsf)+min(tcsf));
        
        [param ffit varacc ffitI domI pk BW] = DoGfit(tcsf,sfdom);
        sfpref(p) = pk;        
        sfmag(p) = (max(ffitI)-min(ffitI))/(max(ffitI)+min(ffitI));
        
        semilogx(sfdom,tcsf,['.' colorid{c} '-']), hold on
        plot(domI,ffitI,'k')
        xlim([min(sfdom) max(sfdom)])
        %ylim([min(ffitI) max(ffitI)])
    end

end


%plot ori curves

figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        %mean over phase
        %kernplot = squeeze(mean(kernC{c}{p},3)); %average over phase
        %kernsigplot = squeeze(mean(kernSigC{c}{p},3)); %average over phase
        kernplot = kernC{c}{p};
        kernsigplot = kernplot;
        
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));  %smooth in all dimensions (used to find slice)
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        
        kernplot = ifftn(fftn(kernplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        
        tcori = squeeze(kernplot(:,bestsfid,tau));
        tcsigori = squeeze(kernsigplot(:,bestsfid,tau));
        
        %tcori = squeeze(mean(kernplot(:,1:2,tau),2));
        %tcsigori = squeeze(mean(kernsigplot(:,1:2,tau),2));
        
        orimuMat{c}(p,:) = tcori;
        orisigMat{c}(p,:) = tcsigori;        
        
        [OMag(p) OAng(p)] = orifind(tcori,oridom);
        
        plot(oridom,tcori,['.' colorid{c} '-']), hold on
    end
    
    %cellS.muBase(p) = mean(cellS.muTimeBlank{p}(tau:tau));
    %cellS.sigBase(p) = mean(cellS.sigTimeBlank{p}(tau:tau));
    
    
end
%%
%Get population dynamics of ori
kernplotmu = 0;
orismoother = getSmoother(1,[.2 1 .2],[1 1],taudom,oridom,sfdom);
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        %mean over phase
        %kernplot = squeeze(mean(kernC{c}{p},3)); %average over phase
        %kernsigplot = squeeze(mean(kernSigC{c}{p},3)); %average over phase
        kernplot = kernC{c}{p};
        %kernsigplot = kernplot;
        
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));  %smooth in all dimensions (used to find slice)
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        
        kernplot = ifftn(fftn(kernplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        %kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [dum bestsfid] = find(ma == max(ma(:)));
%         tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        tau = length(taudom)-2;        
        
        tcori = squeeze(kernplot(:,bestsfid,tau));
        tcoriPre1 = squeeze(kernplot(:,bestsfid,tau-3));
        tcoriPre2 = squeeze(kernplot(:,bestsfid,tau-1));
        tcoriPost = squeeze(kernplot(:,bestsfid,tau+1));  
        zer = squeeze(kernplot(:,bestsfid,end));        
        
        orimuMat{c}(p,:) = tcori/max(tcori);
        orimuMatPre1{c}(p,:) = tcoriPre1/max(tcori);
        orimuMatPre2{c}(p,:) = tcoriPre2/max(tcori);
        orimuMatPost{c}(p,:) = tcoriPost/max(tcori);     
        zermuMat{c}(p,:) = zer/max(tcori);
        
    end
    
    %cellS.muBase(p) = mean(cellS.muTimeBlank{p}(tau:tau));
    %cellS.sigBase(p) = mean(cellS.sigTimeBlank{p}(tau:tau));   
    kernplotmu = kernplotmu + squeeze(kernplot(:,bestsfid,:));
    
end
muPost = fliplr(mean(orimuMatPost{1}));
muma = fliplr(mean(orimuMat{1}));
muPre1 = fliplr(mean(orimuMatPre1{1}));
muPre2 = fliplr(mean(orimuMatPre2{1}));
muZer = fliplr(mean(zermuMat{1}));

stdPost = fliplr(std(orimuMatPost{1}))/sqrt(Ncell);
stdma = fliplr(std(orimuMat{1}))/sqrt(Ncell);
stdPre1 = fliplr(std(orimuMatPre1{1}))/sqrt(Ncell);
stdPre2 = fliplr(std(orimuMatPre2{1}))/sqrt(Ncell);
stdZer = fliplr(std(zermuMat{1}))/sqrt(Ncell);

figure,
errorbar(oridom,muPost,stdPost,'b')
hold on,errorbar(oridom,muma,stdma,'g')
hold on,errorbar(oridom,muPre1,stdPre1,'r')


xlim([-10 180])
xlabel('orientation')

figure,imagesc(oridom, taudom, flipud(kernplotmu'))

%%
%Get population dynamics of sp freq
orismoother = getSmoother(1,[.5 1 .5],1,taudom,oridom,sfdom);
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        %mean over phase
        %kernplot = squeeze(mean(kernC{c}{p},3)); %average over phase
        %kernsigplot = squeeze(mean(kernSigC{c}{p},3)); %average over phase
        kernplot = kernC{c}{p};
        %kernsigplot = kernplot;
        
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));  %smooth in all dimensions (used to find slice)
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother))); %smooth over all dimensions but orientation
        %kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(sfsmoother))); %smooth over all dimensions but orientation
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        %         tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        tau = length(taudom)-2;
        
        tcsf = squeeze(kernplot(bestoriid,:,tau));
        tcsfPre = squeeze(kernplot(bestoriid,:,tau-3));
        tcsfPost = squeeze(kernplot(bestoriid,:,tau+1));        
        
        sfmuMat{c}(p,:) = tcsf/max(tcsf);
        sfmuMatPre{c}(p,:) = tcsfPre/max(tcsf);
        sfmuMatPost{c}(p,:) = tcsfPost/max(tcsf);      
        
    end
    
    %cellS.muBase(p) = mean(cellS.muTimeBlank{p}(tau:tau));
    %cellS.sigBase(p) = mean(cellS.sigTimeBlank{p}(tau:tau));   
    
end
muPost = (mean(sfmuMatPost{1}));
muma = (mean(sfmuMat{1}));
muPre = (mean(sfmuMatPre{1}));
stdPost = (std(sfmuMatPost{1}))/sqrt(Ncell);
stdma = (std(sfmuMat{1}))/sqrt(Ncell);
stdPre = (std(sfmuMatPre{1}))/sqrt(Ncell);

figure,errorbar(sfdom,muPost,stdPost,'b')
hold on,errorbar(sfdom,muma,stdma,'g')
hold on,errorbar(sfdom,muPre,stdPre,'r')
xlabel('spat freq')
%%

cellS.orimuMat = orimuMat{1};
cellS.orisigMat = orisigMat{1};

%Now plot color image of tuning
%OMag = OMag-min(OMag);
OMag = OMag/max(OMag);

mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfmag2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfpref2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

for p = 1:Ncell

    idcell = find(masklabel(:) == celldom(nID(p)));
    
    mag(idcell) = OMag(p);
    ang(idcell) = OAng(p);
    
    sfmag2(idcell) = sfmag(p);
    sfpref2(idcell) = sfpref(p);
    
end

mag = log10(mag+.01);
mag = mag-min(mag(:));
mag = mag/max(mag(:));


figure,
imagesc(ang,'AlphaData',mag,[0 180]), colormap hsv, colorbar, 
hold on
idcirc = [4 7 8 12 13];
for i = 1:length(idcirc)
    p = idcirc(i);
    [idy idx] = find(masklabel == celldom(nID(p)));
    yloc = round(mean(idy)); 
    xloc = round(mean(idx));
    plot(xloc,yloc,'ok','MarkerSize',16)
    hold on
end


figure
imagesc(1:length(sfpref2(1,:)),1:length(sfpref2(:,1)),log10(sfpref2),'Alphadata',sfmag2,log10([sfdom(1) sfdom(end)])), colorbar
for i = 1:length(sfdom)
    domcell{i} = round(sfdom(i)*100)/100;
end
sfvec = round(log10(sfdom)*100)/100;
sfvec(end) = floor(log10(sfdom(end))*100)/100;
sfvec(1) = ceil(log10(sfdom(1))*100)/100;
colorbar('YTick',sfvec,'YTickLabel',domcell)

hold on
idcirc = [4 7 8 12 13];
for i = 1:length(idcirc)
    p = idcirc(i);
    [idy idx] = find(masklabel == celldom(nID(p)));
    yloc = round(mean(idy)); 
    xloc = round(mean(idx));
    plot(xloc,yloc,'ok','MarkerSize',16)
    hold on
end




%%

%time courses
dori = oridom(2)-oridom(1);
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = kernC{c}{p};
    kernplotSig = kernplot;
    
    kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
    kernplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));
    kernplotSig = ifftn(fftn(kernplotSig).*abs(fftn(tsmoother)));
    
    [ma idma] = max(kerndum,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    [mi idmi] = min(kerndum,[],3);
    [worstoriid worstsfid] = find(mi == min(mi(:)));
    
    bestoriid = 4;
    
    orthoriid = bestoriid+round(90/dori);
    if orthoriid>length(oridom)
        orthoriid = orthoriid-length(oridom);
    end

    %x = 7;
    
    tcourse_ma = squeeze(kernplot(bestoriid,bestsfid,:));
    tcourseSig_ma = squeeze(kernplotSig(bestoriid,bestsfid,:));
    tcourse_orth = squeeze(kernplot(orthoriid,bestsfid,:));
    tcourseSig_orth = squeeze(kernplotSig(orthoriid,bestsfid,:));
    tcourse_all = squeeze(mean(mean(kernplot(:,:,:),1),2));
    tcourseSig_all = squeeze(mean(mean(kernplotSig(:,:,:),1),2));

    tcoursemaMat{c}(p,:) = tcourse_ma;
    tcourseorthMat{c}(p,:) = tcourse_orth;
    
    tcourseSigmaMat{c}(p,:) = tcourseSig_ma;
    tcourseSigorthMat{c}(p,:) = tcourseSig_orth;
    
    tcourseAllMat{c}(p,:) = tcourse_all;
    tcourseSigAllMat{c}(p,:) = tcourseSig_all;
    
%     tcourse_ma = mean(kernplot);
%     tcourse_orth = blankresp;
    
    plot(taudom,flipud(tcourse_ma),'.-'), hold on, plot(taudom,flipud(tcourse_orth),'.-r')
    xlim([taudom(1) taudom(end)])

end

mutc = fliplr(mean(tcoursemaMat{1}));
sigtc = fliplr(std(tcoursemaMat{1}))/sqrt(Ncell);
mutc_orth = fliplr(mean(tcourseorthMat{1}));
sigtc_orth = fliplr(std(tcourseorthMat{1}))/sqrt(Ncell);

figure,errorbar(taudom,mutc,sigtc,'k')
hold on
errorbar(taudom,mutc_orth,sigtc_orth,'r')

figure,plot(taudom,flipud(tcoursemaMat{1}'),'k')
hold on, plot(taudom,flipud(tcourseorthMat{1}'),'r')



function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));


function smoother = getSmoother(ktau,kori,ksf,taudom,oridom,sfdom)


ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [kori zeros(1,length(oridom)-length(kori))];
kdum = kori'*ksf;
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));