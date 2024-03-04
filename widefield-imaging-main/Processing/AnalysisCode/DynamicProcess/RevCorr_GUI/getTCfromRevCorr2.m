function getTCfromRevCorr2

%Ian Nauhaus

global cellS DM TC MK

TC = struct;

for i = 1:length(DM.colordom)
    for p = 1:length(cellS.muTime)
        kernC{i}{p} = squeeze(cellS.muTime{p}(:,:,:,i,:));
        kernSigC{i}{p} = squeeze(cellS.sigTime{p}(:,:,:,i,:));
    end    
end

delayWin = [100 500];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-DM.taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-DM.taudom));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions

kernsmooth = getSmoother(100,20,.2,DM.taudom,DM.oridom,DM.sfdom);

%Smoother before taking ori curve
orismoother = getSmoother(10,5,.1,DM.taudom,DM.oridom,DM.sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother(10,10,.1,DM.taudom,DM.oridom,DM.sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother(1,20,.2,DM.taudom,DM.oridom,DM.sfdom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get blank response
ht = zeros(1,length(DM.taudom));
ht(1:3) = [.5 1 .5]; ht = ht/sum(ht);
for c = 1:length(DM.colordom)

    for p = 1:MK.Ncell
        if c == 1
            muBlank{p} = 0;
        end
        if getparam('blankProb')==0
            kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
            kernplot(:,:,:,:) = kernC{c}{p};
            kernplot = mean(kernplot,3); %average over phase
            kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
            muBlankdum = squeeze(mean(kernplot(:,end,:),1)); %average over ori; use last spatial frequency for the blank
        else
            muBlankdum = cellS.muTimeBlank{p};
        end

        muBlankdum = ifft(fft(muBlankdum(:)).*abs(fft(ht(:))));
        muBlank{p} = muBlank{p} + muBlankdum/length(DM.colordom); %average across color
    end
    
end

%subtract the blank
for c = 1:length(DM.colordom)
    for p = 1:MK.Ncell
        for or = 1:length(DM.oridom)
            for sfr = 1:length(DM.sfdom)
                for phs = 1:length(DM.phasedom)
                    kernC{i}{p}(or,sfr,phs,:) = squeeze(kernC{i}{p}(or,sfr,phs,:)) - muBlank{p}(:);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorid = {'r','g','b'};
%get/plot sf curves
%%
figure
for c = 1:length(DM.colordom)
    tcsfall{c} = zeros(MK.Ncell,length(DM.sfdom));
    
    for p = 1:MK.Ncell
        subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
        kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        kernplot = mean(kernplot,3); %average over phase
        kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
      
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));        
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother)));        
        
        [ma idma] = max(kerndum,[],3); %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;       
        
        tcorisfraw = kernplot(:,:,tau);
        
        sfdomI = logspace(log10(DM.sfdom(1)),log10(DM.sfdom(end)-.1),14);
        DM.sfdomI = sfdomI;
        
        tcorisfraw = interp1(DM.sfdom,tcorisfraw',sfdomI,'linear')';
        
%         [param ffitorisf varacc] = Gaussfit2D(1,tcorisfraw,1);
%         [idy idx] = find(ffitorisf == max(ffitorisf(:)));  
%         ffit = ffitorisf(idy,:);
        
        %%%extract sfreq curve using svd
        base = prctile(tcorisfraw(:),20);
        tcorisfraw = phi(tcorisfraw-base);     
        [u s v] = svd(tcorisfraw);                
        tcorisf = u(:,1)*s(1,1)*v(:,1)' + base;   
        [idy idx] = find(tcorisf == max(tcorisf(:)));  %I do it this dumb way because sometimes there are double negatives
        tcsf = tcorisf(idy,:);
        %%%or take slice%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %tcsf = squeeze(kernplot(bestoriid,:,tau));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [param ffit varacc] = DoGfit(tcsf,DM.sfdomI);
          
        TC.tcsfall{c}(p,:) = tcsf;
               
        
        if varacc > .7
            [ma id] = max(ffit);
            TC.sfpref{c}(p) = sfdomI(id);        
            TC.sfmag{c}(p) = (max(ffit)-min(ffit))/(max(ffit)+min(ffit));
            
            [ma idma] = max(ffit); [mi idmi] = min(ffit(idma:end)); idmi = idmi+idma-1;
            thresh = (ma-mi)*.7 + mi;
            [dum fhidum] = min(abs(ffit(idma:idmi) - thresh));
            TC.fhi{c}(p) = sfdomI(fhidum+idma-1);  %high cutoff sfreq
                      
            [dum flodum] = min(abs(ffit(1:idma) - thresh));
            TC.flo{c}(p) = sfdomI(flodum);  %high cutoff sfreq
            
            TC.Qfac{c}(p) = TC.sfpref{c}(p)/(TC.fhi{c}(p)-TC.flo{c}(p));
            TC.sfBW{c}(p) = log2(TC.fhi{c}(p)/TC.flo{c}(p)); %bandwidth in octaves
            
            if TC.sfBW{c}(p) == 0                
                TC.sfBW{c}(p) = NaN;                
            end
            TC.LPness{c}(p) = ffit(1)/ffit(idma);  
        else
            TC.sfpref{c}(p) = NaN;        
            TC.sfmag{c}(p) = 0;
            TC.Qfac{c}(p) = NaN;
            TC.sfBW{c}(p) = NaN;
            TC.LPness{c}(p) = NaN;
            TC.flo{c}(p) = NaN;
            TC.fhi{c}(p) = NaN;
        end
        
%         tcsf = tcsf-min(tcsf);
%         tcsf = tcsf/sum(tcsf);
%         sfpref{c}(p) = sum(tcsf.*DM.sfdom);
        
        semilogx(DM.sfdomI,tcsf,['.' colorid{c} '-']), hold on
        plot(DM.sfdomI,ffit,'k')       
        axis tight
% 
%         title(num2str(TC.Qfac{c}(p)))     
        
    end

end
if ~isempty(squeeze(find(TC.Qfac{1}>3 & TC.Qfac{1}<inf)))
    'hello'
end
%%
%get/plot ori curves
colorvec = zeros(MK.Ncell,length(DM.colordom));
figure
for c = 1:length(DM.colordom)
    
    tcoriall{c} = zeros(MK.Ncell,length(DM.oridom));
    for p = 1:MK.Ncell
        subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
    
        %mean over phase
        kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        kernplot = mean(kernplot,3); %average over phase
        kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
        
        kernsigplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernsigplot(:,:,:,:) = kernSigC{c}{p};
        kernsigplot = mean(kernsigplot,3); %average over phase
        kernsigplot = reshape(kernsigplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
        
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));  %smooth in all dimensions (used to find slice)
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window        
        kernplot = ifftn(fftn(kernplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(orismoother))); %smooth over all dimensions but orientation
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestsfid) + delayWinID(1) - 1;
        
        %tcori = squeeze(kernplot(:,bestsfid,tau));
        %tcsigori = squeeze(kernsigplot(:,bestsfid,tau));
        
        tcorisfraw = kernplot(:,:,tau);
        
        dI = 5; DM.oridomI = 0:dI:180-dI;
        tcorisfraw = interp1([DM.oridom 180], [tcorisfraw; tcorisfraw(1,:)], [DM.oridomI 180]);
        tcorisfraw = tcorisfraw(1:end-1,:);
        
%         [param ffitorisf varacc] = Gaussfit2D(1,tcorisfraw,1);
%         param(1:4) = param(1:4)*dI;
%         param(1) = param(1)+DM.oridom(1);
%         
%         [idy idx] = find(ffitorisf == max(ffitorisf(:)));  
%         ffit = ffitorisf(:,idx);
        
        %%%extract ori curve using svd
        base = prctile(tcorisfraw(:),20);
        tcorisfraw = phi(tcorisfraw-base);
        [u s v] = svd(tcorisfraw);
        tcori = u(:,1)*s(1,1)*v(:,1)' + base;
        [idy idx] = find(tcori == max(tcori(:)));  %I do it this dumb way because sometimes there are double negatives
        tcori = tcori(:,idx);
        %%%or take slice %%%%%%%%%%%%%%%%
        %tcori = squeeze(mean(kernplot(:,1:2,tau),2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        [param ffit varacc] = Gaussfit(DM.oridomI,tcori',1);
   
        tcsigori = squeeze(mean(kernsigplot(:,1:2,tau),2));         
             
        [TC.OMag{c}(p) TC.OAng{c}(p)] = orifind(tcori,DM.oridomI); 
        
        if varacc < .7
            param = param*NaN;
            OMag{c}(p) = 0;
        end
        TC.opref{c}(p) = param(1);
        TC.orisig{c}(p) = param(2);       
        
        plot(DM.oridomI,tcori,['.' colorid{c} '-']), hold on
          
        plot(DM.oridomI,ffit,'k')
        title([num2str(round(TC.OMag{c}(p)*100)/100) '  ' num2str(TC.orisig{c}(p))])

        colorvec(p,c) = param(3);
        
        TC.tcoriall{c}(p,:) = tcori;
        orimuMat{c}(p,:) = tcori;
        orisigMat{c}(p,:) = tcsigori;  
        
    end
    
    if ~isempty(cellS.muTimeBlank{p})
        cellS.muBase(p) = mean(cellS.muTimeBlank{p}(tau:tau));
        cellS.sigBase(p) = mean(cellS.sigTimeBlank{p}(tau:tau));
    end
    
    TC.OMag{c} = phi(TC.OMag{c});
end

cellS.orimuMat = orimuMat{1};
cellS.orisigMat = orisigMat{1};


%control
% [dum id] = sort(rand(size(TC.OAng{1})));
% TC.OAng{1} = TC.OAng{1}(id);

%% Now phase

figure
for p = 1:MK.Ncell
    subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)

    for c = 1:length(DM.colordom)
        
        kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        
        kerndum = kernplot;
        %kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));        
        tcourse = squeeze(mean(mean(mean(kerndum,1),2),3));
        [dum id] = max(tcourse);
        kerndum = kernplot(:,:,:,id);
        kerndum = reshape(kerndum,length(DM.oridom),length(DM.sfdom),length(DM.phasedom)); %squeeze time dimension
        
        [ma idma] = max(kerndum,[],3);  %find maxima from smoothed version
        [bestoriid bestsfid] = find(ma == max(ma(:)));
        
        tcphase = squeeze(kerndum(bestoriid,bestsfid,:));
        
        f1 = abs(sum(tcphase.*exp(1i*DM.phasedom'*pi/180)));
        f0 = mean(tcphase);
        
        TC.F1F0{c}(p) = f1/f0;
        
        %The Nishimoto et al version
%         f1 = abs(tcphase(1)-tcphase(3)) + abs(tcphase(2)-tcphase(4));
%         f0 = sum(tcphase);
%         TC.F1F0{c}(p) = 2*f1/f0;
        
        plot(DM.phasedom,tcphase,['.' colorid{c} '-']), hold on
        
        title(num2str(TC.F1F0{c}(p)))
        
    end

end

id = find(TC.F1F0{c}<0);
TC.F1F0{c}(id) = 0;
id = find(TC.F1F0{c}>2);
TC.F1F0{c}(id) = 2;
figure,hist(TC.F1F0{c},linspace(0,2,20))


function smoother = getSmoother(tausig,orisig,sfsig,taudom,oridom,sfdom)

if length(sfdom) == 1
    ksf = 1;
end
if length(oridom) == 1
    kori = 1;
end

ktau = exp(-(taudom-mean(taudom)).^2/(2*tausig^2)); %tausig is in ms

dom = log2(sfdom)-mean(log2(sfdom));
ksf = exp(-dom.^2/(2*sfsig^2)); %sfsig is in octaves

oridomdum = linspace(0,180,length(oridom)+1);
oridomdum = oridomdum(1:end-1);  %I use this one in case oridom wraps around
kori = exp(-(oridomdum-oridomdum(ceil(end/2))).^2/(2*orisig^2)); %orisig is in degrees

kdum = kori'*ksf;
smoother = zeros(length(kori),length(ksf),length(ktau));
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));


function tcori = intCirc(oridom,tcori,oridomI)

Ifac = length(oridomI)/length(oridom);

[dum idma] = max(tcori);
shifta = round(length(tcori)/2)-idma;
tcori = circshift(tcori(:)', [0 shifta]);

tcori = interp1(oridom,tcori,oridomI,'spline');

shiftb = round(-shifta*Ifac);
tcori = circshift(tcori, [0 shiftb]);
