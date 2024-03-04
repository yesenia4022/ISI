function Gkernelplots2

global ACQinfo Analyzer cellS maskS G_RChandles

%%%%

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
[nID] = getNeuronMask;

%%%%

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string') ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with dtau spacing, and end at an estimate of kernDel(2)

delayWin = [250 600];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

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
    for p = 1:length(cellS.muTime)
        kernC{i}{p} = squeeze(cellS.muTime{p}(:,:,:,i,:));
        kernSigC{i}{p} = squeeze(cellS.sigTime{p}(:,:,:,i,:));
    end    
end

Ncell = length(cellS.muTime);
NT = getnotrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)

kernsmooth = getSmoother([.1 .5 1 .5 .1],[.5 1 .5],[.1 1 .1],taudom,oridom,sfdom);

%Smoother before taking ori curve
orismoother = getSmoother([.1 .5 1 .5 .1],1,[.2 1 .2],taudom,oridom,sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother([.2 .5 1 .5 .2],[.3 1 .3],1,taudom,oridom,sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother(1,[.3 1 .3],[.1 1 .1],taudom,oridom,sfdom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorid = {'r','g','b'};
%plot sf curves

figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        kernplot = kernC{c}{p};
        kernplot = squeeze(mean(kernplot,3)); %average over phase
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother)));        
        
        [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        
        tcsf = squeeze(kernplot(bestoriid,:,tau));
        
%         [dum sfpref(p)] = max(tcsf);
%         sfpref(p) = sfdom(sfpref(p));
%         sfmag(p) = (max(tcsf)-min(tcsf));
        
        [param ffit varacc ffitI domI pk BW] = DoGfit(tcsf,sfdom);
        sfpref(p) = pk;        
        sfmag(p) = (max(ffitI)-min(ffitI))/(max(ffitI)+min(ffitI));
        
        semilogx(sfdom,tcsf,['.' colorid{c} '-']), hold on
        plot(domI,ffitI,'k')
        xlim([min(sfdom) max(sfdom)])
        ylim([min(ffitI) max(ffitI)])
    end

end


%plot ori curves

figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        %mean over phase
        kernplot = squeeze(mean(kernC{c}{p},3)); %average over phase
        kernsigplot = squeeze(mean(kernSigC{c}{p},3)); %average over phase
        
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
    
    cellS.muBase(p) = mean(cellS.muTimeBlank{p}(tau:tau));
    cellS.sigBase(p) = mean(cellS.sigTimeBlank{p}(tau:tau));
    
    
end

cellS.orimuMat = orimuMat{1};
cellS.orisigMat = orisigMat{1};


figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    kerndum = kernplot;
    %kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    tcourse = squeeze(mean(mean(mean(kerndum,1),2),3));
    [dum id] = max(tcourse);
    kerndum = squeeze(kernplot(:,:,:,id));
    
    [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    for c = 1:length(colordom)
        
        kernplot = kernC{c}{p};
        kerndum = kernplot;
        %kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));        
        tcourse = squeeze(mean(mean(mean(kerndum,1),2),3));
        [dum id] = max(tcourse);
        kerndum = squeeze(kernplot(:,:,:,id));
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        
        tcphase = squeeze(kerndum(bestoriid,bestsfid,:));
        
        plot(phasedom,tcphase,['.' colorid{c} '-']), hold on
    end

end

figure,scatter(OAng,sfpref,'.'), xlabel('pref orientation'), ylabel('pref spatial frequency')


%Now plot color image of tuning
OMag = OMag-min(OMag);
OMag = OMag/max(OMag);

mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfmag2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfpref2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
for p = 2:Ncell

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
imagesc(ang,'AlphaData',mag,[0 180]), colormap hsv, colorbar
figure
imagesc(log10(sfpref2),'Alphadata',sfmag2,log10([sfdom(1) sfdom(end)])), colorbar


%time courses
figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = cellS.muTime{p};
    kernplotSig = cellS.sigTime{p};
    kernplot = squeeze(mean(kernplot,3)); %average over phase    
    kernplotSig = squeeze(mean(kernplotSig,3));
    kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
    kernplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));
    kernplotSig = ifftn(fftn(kernplotSig).*abs(fftn(tsmoother)));
    
    [ma idma] = max(kerndum,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    [mi idmi] = min(kerndum,[],3);
    [worstoriid worstsfid] = find(mi == min(mi(:)));

    %x = 7;

    tcourse_ma = squeeze(kernplot(bestoriid,bestsfid,:));
    tcourseSig_ma = squeeze(kernplotSig(bestoriid,bestsfid,:));
    tcourse_orth = squeeze(kernplot(worstoriid,worstsfid,:));
    tcourseSig_orth = squeeze(kernplotSig(worstoriid,worstsfid,:));
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
    
    plot(taudom,tcourse_ma,'.-'), hold on, plot(taudom,tcourse_orth,'.-r')
    xlim([taudom(1) taudom(end)])

end

cellS.tcoursemaMat = tcoursemaMat{1};
cellS.tcourseorthMat = tcourseorthMat{1};
cellS.tcourseSigmaMat = tcourseSigmaMat{1};
cellS.tcourseSigorthMat = tcourseSigorthMat{1};
cellS.tcourseAllMat = tcourseAllMat{1};
cellS.tcourseSigAllMat = tcourseSigAllMat{1};




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