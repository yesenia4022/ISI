function Gplotlogmap_RC

global ACQinfo Analyzer cellS maskS G_RChandles

%%%%

tauN = str2num(get(G_RChandles.kernelLength,'string'));

masklabel = bwlabel(maskS.bwCell{1},4);
celldom = unique(masklabel);
nID = getNeuronMask;

%%%%

eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)

delayWin = [100 1000];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
%load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
load(['I:\neurostuff\log_files\' Analyzer.M.anim '\' expt])

[domains seqs] = getSeqInfo;


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

Ntau = Ntau+1;

%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)

kernsmooth = getSmoother([.1 .5 1 .5 .1],[.1 1 .1],[.5 1 .5],taudom,oridom,sfdom);

%Smoother before taking ori curve
orismoother = getSmoother([.1 .5 1 .5 .1],[.1 1 .1],1,taudom,oridom,sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother([.2 .5 1 .5 .2],1,[.3 1 .3],taudom,oridom,sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother(1,[.4 1 .4],[.1 1 .1],taudom,oridom,sfdom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorid = {'r','g','b'};
%plot sf curves

%figure
for p = 1:Ncell
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    for c = 1:length(colordom)
        kernplot = kernC{c}{p};
        kernplot = squeeze(mean(kernplot,3)); %average over phase
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
        kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
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
        ylim([min(tcsf) max(tcsf)])
    end

end


%time courses
%figure
for p = 1:Ncell
    %subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    kernplot = cellS.muTime{p};
    kernplotSig = cellS.sigTime{p};
    kernplot = squeeze(mean(kernplot,3)); %average over phase    
    kernplotSig = squeeze(mean(kernplotSig,3));
    kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
    kernplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));
    kernplotSig = ifftn(fftn(kernplotSig).*abs(fftn(tsmoother)));
    
    %Reduce error bars based on the smoothing operation
    tsmoothSlice = tsmoother(:,:,1);
    D = sum(tsmoothSlice(:)/max(tsmoothSlice(:)));
    kernplotSig = kernplotSig/sqrt(D);
    
    [ma idma] = max(kerndum,[],3);
    [bestoriid bestsfid] = find(ma == max(ma(:)));
    [mi idmi] = min(kerndum,[],3);
    [worstoriid worstsfid] = find(mi == min(mi(:)));

    %x = 7;

    tcourse_ma = squeeze(kernplot(bestoriid,bestsfid,:));
    tcourseSig_ma = squeeze(kernplotSig(bestoriid,bestsfid,:));
    tcourse_orth = squeeze(kernplot(worstoriid,worstsfid,:));
    tcourseSig_orth = squeeze(kernplotSig(worstoriid,worstsfid,:));

    tcoursemaMat{c}(p,:) = tcourse_ma;
    tcourseorthMat{c}(p,:) = tcourse_orth;
    
    tcourseSigmaMat{c}(p,:) = tcourseSig_ma;
    tcourseSigorthMat{c}(p,:) = tcourseSig_orth;
    
%     tcourse_ma = mean(kernplot);
%     tcourse_orth = blankresp;
    
    %plot(taudom,tcourse_ma), hold on, plot(taudom,tcourse_orth,'-r')
    %xlim([taudom(1) taudom(end)])

end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


%Now plot color image of tuning

sfmag2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
sfpref2 = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
for p = 2:Ncell

    idcell{p} = find(masklabel(:) == celldom(nID(p)));
    
    sfmag2(idcell{p}) = sfmag(p);
    sfpref2(idcell{p}) = sfpref(p);
    
end

sfmag2 = sfmag2(3:end-2,3:end-2); sfpref2 = sfpref2(3:end-2,3:end-2);

sfmag2 = sfmag2-min(sfmag2(:));
sfmag2 = sfmag2/max(sfmag2(:));

dim = size(sfmag2);
set(gcf,'Color',[1 1 1]);
%%
anatflag = 1;
figure
if anatflag
    
    imanat = maskS.im{1}(3:end-2,3:end-2);
    
    mi = prctile(imanat(:),1);    
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),99);
    imanat = imanat/ma;

    imfunc = sfpref2;
    imfunc = imfunc-sfdom(1);
    imfunc = imfunc/(sfdom(end)-sfdom(1));
    imfunc = round(imfunc*63+1);
    
    id = find(imfunc<=0);
    imfunc(id) = 1;

    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = sfmag2(i,j)*jetid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+sqrt(imanat);

    imout = imout/max(imout(:));
    
    x = image(imout,'CDataMapping','direct');

else
    
    imfunc = sfpref2;
    imfunc = imfunc-log10(logdom(1));
    imfunc = imfunc/(log10(logdom(end))-log10(logdom(1)));
    imfunc = round(imfunc*63+1);
    
    image(1:length(pref(1,:)),1:length(pref(:,1)),imfunc,CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    colorbar
end


for i = 1:length(sfdom)
    domcell{i} = round(sfdom(i)*10)/10;
end
iddom = linspace(1,64,length(sfdom));
colorbar('YTick',iddom,'YTickLabel',domcell)

axis image

fh = gcf;
%%
%%%%%%%%%%%%%%%

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

global cellS ACQinfo maskS G_RChandles

pos = round(get(event_obj,'Position')); %pos(1) is column dimension

masklabel = bwlabel(maskS.neuronmask(3:end-2,3:end-2),4);

cellID = masklabel(pos(2),pos(1));

oridom = linspace(0,180,length(cellS.orimuMat(1,:))+1);
oridom = oridom(1:end-1);

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string') ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod);
taudom = linspace(kernDel(1),kernDel(2),Ntau);

if cellID ~= 0 %if not neuropil

    figure(98)
    subplot(1,2,1)

    tc = cellS.tcoursemaMat(cellID,:);
    tcSig = cellS.tcourseSigmaMat(cellID,:);
    fill([taudom fliplr(taudom)],[tc-tcSig fliplr(tc+tcSig)],'r')
    hold on
    plot(taudom,tc,'k','LineWidth',2)
    xlim([taudom(1) taudom(end)])
    ylabel('Z')
    xlabel('ms')
    title(['Location x' num2str(pos(1)) ' y' num2str(pos(2))])
    hold off
    
    subplot(1,2,2)
    
    tc = cellS.orimuMat(cellID,:);
    [dum ind] = max(tc);
    tc = circshift(tc,[0 round(length(tc)/2)-ind+1]);
    
    domI = linspace(oridom(1),oridom(end),4*length(oridom));
    [tcI] = interp1(oridom,tc,domI,'spline');

    [param ffit varacc ffitI domIfit] = Gaussfit(domI,tcI,0);
    
    %errorbar(oridom,cellS.orimuMat(cellID,:),cellS.orisigMat(cellID,:))  
    
    plot(domIfit,ffitI,'k')
    hold on
    plot(oridom,tc,'.r')
    title(num2str(param(2)))
    hold off
end


tar = get(get(event_obj,'Target'));
data = tar.CData;

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Ori: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180) ' deg']};

  
   
function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));


function smoother = getSmoother(ktau,ksf,kori,taudom,oridom,sfdom)

ktau = [ktau zeros(1,length(taudom)-length(ktau))];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [kori zeros(1,length(oridom)-length(kori))];
kdum = kori'*ksf;
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));

