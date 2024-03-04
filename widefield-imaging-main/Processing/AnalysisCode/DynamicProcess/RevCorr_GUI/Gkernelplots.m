function [oriPop] = Gkernelplots(tauN)

global ACQinfo Analyzer maskS cellS

%%%%

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
dtau = acqPeriod;
Ntau = round(tauN/acqPeriod);
taudom = 0:dtau:dtau*Ntau;

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])

%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;
colordom = domains.colordom;

for i = 1:length(colordom)
    for p = 1:length(cellS.mu)
        kernC{i}{p} = squeeze(cellS.mu{p}(:,:,:,i,:));
        kernSigC{i}{p} = squeeze(cellS.sig{p}(:,:,:,i,:));
    end    
end

Ncell = length(cellS.mu);
NT = getnotrials;

Ntau = Ntau+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)
ktau = [.1 .5 1 .5 .1];
ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = [.2 1 .2];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [.5 1 .5];
kori = [kori zeros(1,length(oridom)-length(kori))];

kdum = kori'*ksf;
for i = 1:length(ktau)
    kernsmooth(:,:,i) = kdum*ktau(i);
end
kernsmooth = kernsmooth/sum(kernsmooth(:));

%Smoother before taking ori curve
ktau = [.1 .5 1 .5 .1];
ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = [.2 1 .2];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [1];
kori = [kori zeros(1,length(oridom)-length(kori))];
kdum = kori'*ksf;
for i = 1:length(ktau)
    orismoother(:,:,i) = kdum*ktau(i);
end
orismoother = orismoother/sum(orismoother(:));

%Smoother before taking sf curve
ktau = [.2 .5 1 .5 .2];
ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = 1;
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [.3 1 .3];
kori = [kori zeros(1,length(oridom)-length(kori))];
kdum = kori'*ksf;
for i = 1:length(ktau)
    sfsmoother(:,:,i) = kdum*ktau(i);
end
sfsmoother = sfsmoother/sum(sfsmoother(:));

%Smoother before taking temporal curve
ktau = 1;
ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = [.2 1 .2];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [.3 1 .3];
kori = [kori zeros(1,length(oridom)-length(kori))];
kdum = kori'*ksf;
for i = 1:length(ktau)
    tsmoother(:,:,i) = kdum*ktau(i);
end
tsmoother = tsmoother/sum(tsmoother(:));

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
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother)));
        
        [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid);
        
        tcsf = squeeze(kernplot(bestoriid,:,tau));
        
        [dum sfpref(p)] = max(tcsf);
        sfpref(p) = sfdom(sfpref(p));
        sfmag(p) = (max(tcsf)-min(tcsf));
        
        semilogx(sfdom,tcsf,['.' colorid{c} '-']), hold on
        xlim([min(sfdom) max(sfdom)])
        %ylim([min(tcsf) max(tcsf)])
    end

end


%plot ori curves

figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        kernplot = kernC{c}{p};
        kernplot = squeeze(mean(kernplot,3)); %average over phase
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));  %smooth in time
        kernplot = ifftn(fftn(kernplot).*abs(fftn(orismoother))); %smooth over orientation
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid);
        
        %tcori = squeeze(kernplot(:,bestsfid,tau));
        
        tcori = squeeze(mean(kernplot(:,1:2,tau),2));
        orimuPop{c}(p,:) = tcori;
        
        [OMag(p) OAng(p)] = orifind(tcori,oridom);
        
        plot(oridom,tcori,['.' colorid{c} '-']), hold on
    end

end
oriPop = oriPop{c};

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


%Now plot color image of tuning
OMag = OMag-min(OMag);
OMag = OMag/max(OMag);

mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfmag2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfpref2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
for p = 2:Ncell
    
    idcell = find(masklabel(:) == celldom(p));
    
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

    kernplot = cellS.mu{p};
    kernplot = squeeze(mean(kernplot,3)); %average over phase
    blankresp = cellS.mu{p}(end,:);
    
    kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    kernplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));
    
    [ma idma] = max(kerndum,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    [mi idmi] = min(kerndum,[],3);
    [worstoriid worstsfid] = find(mi == min(mi(:)));

    %x = 7;

    tcourse_ma = squeeze(kernplot(bestoriid,bestsfid,:));   
    tcourse_orth = squeeze(kernplot(worstoriid,worstsfid,:)); 
    
%     tcourse_ma = mean(kernplot);
%     tcourse_orth = blankresp;
    
    plot(taudom,tcourse_ma,'.-'), hold on, plot(taudom,tcourse_orth,'.-r')
    xlim([taudom(1) taudom(end)])

end





function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));