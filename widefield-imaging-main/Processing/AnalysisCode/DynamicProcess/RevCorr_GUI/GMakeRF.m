function RF = GMakeRF(kern,tauN)

global Analyzer ACQinfo G_RChandles DM

%kern is generated from getTCfromRevCorr4

%%%%

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt],'frate')

if ~exist('frate')
    frate = 60.0044;
end

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file)

trialdom = 1:getnotrials;
[domains seqs] = getSeqInfo(trialdom);


%%%%%%%%%%%%%%%%%%%

oridom = domains{trialdom(1)}.oridom;
sfdom = domains{trialdom(1)}.sfdom;
phasedom = domains{trialdom(1)}.phasedom;
colordom = domains{trialdom(1)}.colordom;

Ncell = length(kern{1});
NT = getnotrials;

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
dtau = acqPeriod;
Ntau = round(tauN/acqPeriod);

taudom = DM.taudom;
delayWin = [150 350];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-DM.taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-DM.taudom));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot sf curves

stimsize = getparam('x_size');

Wcm = stimsize/360*(2*pi*Analyzer.M.screenDist);
pixpercm = Analyzer.M.xpixels/Analyzer.M.screenXcm;
Npix = Wcm*pixpercm;

stimsize = 1/(sfdom(1));
%stimsize = getparam('x_size')
%phi = linspace(-stimsize/2,stimsize/2,Npix);
phi = linspace(0,stimsize,Npix);
colors = {'r','g','b'};
%%
Dec = 1;
RFall = 0;
clear RFiradon
figure
for p = 1:Ncell
    
    for c = 1:length(colordom)
        kernplot = kern{c}{p};
        
        %kernplot = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
        
        Frespdum = squeeze(mean(kernplot(:,:,:,delayWinID(1):delayWinID(2)),4));
        %Frespdum = Frespdum - min(Frespdum(:));
        %Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2)-Frespdum(:,:,3)-Frespdum(:,:,4));
        %Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2));
        %Frespdum = randn(size(Frespdum));
        
        
        Frespdum = Frespdum-prctile(Frespdum(:),50);
        id = find(Frespdum(:) < 0);
        Frespdum(id) = 0;
        
      
        
        for ori = 1:length(oridom)
            
            RF = 0;
            for sf = 1:length(sfdom)
                
                for phase = 1:4
                    
                    xphi = sfdom(sf)*phi*2*pi/Dec; %radians
                    
                    RF = RF + Frespdum(ori,sf,phase)*cos(xphi-phasedom(phase)*pi/180);
          
                end
                
            end
            RFiradon(:,ori) = RF;  %projection for each exis
        end
        
        %tc = mean(Frespdum,3);
%         tc = squeeze(mean(Frespdum(:,2:4),2));
        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
        
        %[dum id] = max(var(RFiradon));
        %plot(RFiradon(:,id))
        %ylim([-2 2])
        
        RF = iradon(RFiradon,oridom);
        %hold on
        imagesc(RFiradon)
        
        %plot(mean(RFiradon,2),colors{c})
        
        %plot(squeeze(Frespdum(2,2,:)))
        
        RFall = RFall+RF;
        
        
        imagesc(RF,[prctile(RF(:),1) prctile(RF(:),99)]), axis image
        colorbar
        
        %axis off
        drawnow
    end
    
end

%figure,imagesc(RFall)