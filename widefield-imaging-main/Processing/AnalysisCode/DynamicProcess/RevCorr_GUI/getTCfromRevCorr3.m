function getTCfromRevCorr3

%Ian Nauhaus

%3 is from 2, and it combines the spatial frequency and orientation loops

global cellS DM TC MK idExamp kernC kernSigC

TC = struct;

kernC = cell(1,length(DM.colordom));

for i = 1:length(DM.colordom)
    for p = 1:length(cellS.muTime)
        kernC{i}{p} = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernSigC{i}{p} = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        
        %Need to preserve dimension, even if it is only one element (i.e. one orientation)
        kernC{i}{p}(:,:,:,:) = squeeze(cellS.muTime{p}(:,:,:,i,:));
        kernSigC{i}{p}(:,:,:,:) = squeeze(cellS.sigTime{p}(:,:,:,i,:));
    end    
end

%t2 = getparam('h_per')*10+150; %Keep it restricted to the "onset response" ??
t2 = 500;

delayWin = [100 t2];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-DM.taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-DM.taudom));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions

tdom = DM.taudom(1:end);

kernsmooth = getSmoother(30,10,.2,tdom,DM.oridom,DM.sfdom); %for establising time-to-peak

%Smoother before taking ori curve
%orismoother = getSmoother(10,5,.1,DM.taudom,DM.oridom,DM.sfdom);

%Smoother before taking sf curve
sfsmoother = getSmoother(50,5,.05,tdom,DM.oridom,DM.sfdom);

%Smoother before taking temporal curve
tsmoother = getSmoother(30,10,.05,tdom,DM.oridom,DM.sfdom);  %this is used to determine a "significant" response, 
                                                                %so keep smoothing minimal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get blank response
ht = zeros(1,length(DM.taudom));
ht(1:3) = [.3 1 .3]; ht = ht/sum(ht);
for c = 1:length(DM.colordom)

    for p = 1:MK.Ncell
        if c == 1
            muBlank{p} = 0;
            sigBlank{p} = 0;
        end
       if getparam('blankProb')==0
            kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
            kernplot(:,:,:,:) = kernC{c}{p};
            kernplot = mean(kernplot,3); %average over phase
            kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
            muBlankdum = squeeze(mean(kernplot(:,end,:),1)); %average over ori; use last spatial frequency for the blank
            
            kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
            kernplot(:,:,:,:) = kernSigC{c}{p};
            kernplot = mean(kernplot,3); %average over phase
            kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
            sigBlankdum = squeeze(mean(kernplot(:,end,:),1)); %average over ori; use last spatial frequency for the blank
        else
            muBlankdum = cellS.muTimeBlank{p};
            sigBlankdum = cellS.sigTimeBlank{p};
        end

        muBlankdum = ifft(fft(muBlankdum(:)).*abs(fft(ht(:))));
        muBlank{p} = muBlank{p} + muBlankdum/length(DM.colordom); %average across color
        
        sigBlankdum = ifft(fft(sigBlankdum(:)).*abs(fft(ht(:))));
        sigBlank{p} = sigBlank{p} + sigBlankdum/length(DM.colordom); %average across color
    end
    
end

%subtract the blank
for c = 1:length(DM.colordom)
    for p = 1:MK.Ncell
        
        [dum idz] = min(abs(DM.taudom-0));
        
        base = kernC{c}{p}(:,:,:,idz);
        base = mean(base(:));
        
        muBlankdum = muBlank{p}(:) - muBlank{p}(idz) + base; %give them the same value at t = 0;
        for or = 1:length(DM.oridom)
            for sfr = 1:length(DM.sfdom)
                for phs = 1:length(DM.phasedom)
                    kernC{c}{p}(or,sfr,phs,:) = squeeze(kernC{c}{p}(or,sfr,phs,:)) - muBlankdum;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorid = {'r','g','b'};
%get sf and ori curves
%%
figure
for c = 1:length(DM.colordom)
    idex = 1;
    for p = 1:MK.Ncell
        subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
        kernplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernplot(:,:,:,:) = kernC{c}{p};
        kernplot = mean(kernplot,3); %average over phase
        kernplot = reshape(kernplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
        
        kernsigplot = zeros(length(DM.oridom),length(DM.sfdom),length(DM.phasedom),length(DM.taudom));
        kernsigplot(:,:,:,:) = kernSigC{c}{p};
        kernsigplot = mean(kernsigplot,3); %average over phase
        kernsigplot = reshape(kernsigplot,length(DM.oridom),length(DM.sfdom),length(DM.taudom)); %squeeze phase dimension
        
        %kernplot = diff(kernplot,[],3);
        %kernsigplot = diff(kernsigplot,[],3);
 
        %compute optimal time delay
        kerndum = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));        
        kerndum = kerndum(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window   
        maxprof = max(reshape(kerndum(:),length(DM.oridom)*length(DM.sfdom),length(delayWinID(1):delayWinID(2))));
        [dum idma] = max(maxprof); %time delay with max variance
        TC.tauID{c}(p) = idma + delayWinID(1) - 1;
        
        %new sigma after smoothing requires the following
        %Etc2 = kernsigplot.^2 + kernplot.^2;       %get E(x^2) = var(x) + (E(x))^2 
        %kernsigplot = ifftn(fftn(Etc2.^2).*abs(fftn(sfsmoother.^2)));
        %kernsigplot = kernsigplot - (ifftn(fftn(kernplot).*abs(fftn(sfsmoother)))).^2; %This is what my math said
        %kernsigplot = sqrt(kernsigplot);
        %kernsigplot = ifftn(fftn(kernsigplot).*abs(fftn(sfsmoother)));
        
        %These are for plotting the time courses in the examples and for
        %establishing significance
        kernTplot = ifftn(fftn(kernplot).*abs(fftn(tsmoother)));  
        kernsigTplot = ifftn(fftn(kernsigplot).*abs(fftn(tsmoother)));
        
        %This is for making the tuning curves
        kernplot = ifftn(fftn(kernplot).*abs(fftn(sfsmoother))); %make sure this is done after computing kernsigplot
             
        tcorisfraw = kernplot(:,:,TC.tauID{c}(p));        
        tcorisfrawsig = kernsigplot(:,:,TC.tauID{c}(p)); 
        
        w = phi(mean(tcorisfraw,2))+eps;  %estimate of ori curve 
        %[dum w] = Gaussfit(DM.oridom,w',1); w = w';
        %w = phi(w-max(w(:))/2);
        w = w/sum(w(:));
        tcsf = w'*tcorisfraw; %weighted average over ori
        %tcsf = phi(tcsf);
        
        w = phi(mean(tcorisfraw,1))+eps;  %estimate of spat freq
        %[dum w] = DoGfit(w,DM.sfdom);
        %w = phi(w-max(w(:))/2);
        w = w/sum(w(:));
        tcori = tcorisfraw*w';
        %tcori = phi(tcori);
        
        DM.sfdomI = logspace(log10(DM.sfdom(1)),log10(DM.sfdom(end)-.1),11);     %first interpolate in sp freq   
        %DM.sfdomI = linspace(DM.sfdom(1),DM.sfdom(end)-.1,14); 
        tcsfI = interp1(DM.sfdom',tcsf',DM.sfdomI','spline')';        

        
        %[RFenv RFsize] = getRFenvelope(tcsf,DM.sfdom); 
        RFenv = NaN; RFsize = NaN;
        
        if length(DM.oridom)>1
            dI = 5; DM.oridomI = 0:dI:180-dI;
            tcoriI = interp1([DM.oridom 180], [tcori; tcori(1)], [DM.oridomI 180]);  %interpolate in ori
            tcoriI = tcoriI(1:end-1);
        else
            DM.oridomI = DM.oridom;
        end        
        
        if isnan(norm(tcoriI)) || isnan(norm(tcsfI))
            'hello'
        end
        
        %[paramdum ffit_sf{p} varacc_sf ffitIsf domIsf] = DoGfit(tcsfI,DM.sfdomI);
        [paramdum ffit_sf{p} varacc_sf ffitIsf domIsf] = SoGfit(DM.sfdomI,tcsfI); sfBWLin = 2*paramdum(2);
        
        %[RFenv RFsize] = getRFenvelope(tcsf,DM.sfdom); 
        
        if length(DM.oridomI)>1
            [param ffit_ori{p} varacc_ori sigma] = Gaussfit(DM.oridomI,tcoriI,1);  param(2) = sigma;
            %[mu sigma ffit_ori{p} varacc_ori] = CircGaussFit(tcoriI);    param(1) = mu; param(2) = sigma;
        else
            param = NaN*ones(1,4);
            ffit_ori{p} = NaN;
            varacc_ori = 1;
        end

        [idori idsf] = find(tcorisfraw == max(tcorisfraw(:)));
        tcoursema = squeeze(kernTplot(idori,idsf,:));
        tsigcoursema = squeeze(kernsigTplot(idori,idsf,:));

        [idori idsf] = find(tcorisfraw == min(tcorisfraw(:)));
        tcoursemi = squeeze(kernplot(idori,idsf,:));
        tsigcoursemi = squeeze(kernsigplot(idori,idsf,:));
        
        [dum idma] = max(tcoursema);
        %TC.SNR{c}(p) = (tcoursema(idma) - tcoursema(idz))/(tsigcoursema(idma) + tsigcoursema(idz));
        %TC.SNR{c}(p) = (tcoursema(idma))/(tsigcoursema(idma) + sigBlank{p}(idma));
        TC.SNR{c}(p) = (tcoursema(idma) - tcoursemi(idma))/(tsigcoursema(idma) + tsigcoursemi(idma));
        
        TC.tcsfall{c}(p,:) = tcsfI;
        TC.tcoriall{c}(p,:) = tcoriI;
        TC.tcsfall_fit{c}(p,:) = ffit_sf{p};
        TC.tcoriall_fit{c}(p,:) = ffit_ori{p};               
        
        if varacc_ori > .6 & varacc_sf > .6 & TC.SNR{c}(p) > 1
            
            %First get all the spatial frequency metrics
            [ma id] = max(ffitIsf);
            TC.sfpref{c}(p) = domIsf(id);        
            TC.sfmag{c}(p) = (max(ffitIsf)-min(ffitIsf))/(max(ffitIsf)+min(ffitIsf));
            TC.sfmag{c}(p) = phi(TC.sfmag{c}(p));  %very rarely is it <0
            
            [ma idma] = max(ffitIsf); [mi idmi] = min(ffitIsf(idma:end)); idmi = idmi+idma-1;  %better to use high sfreq as baseline
            thresh = (ma-mi)/sqrt(2) + mi;
            [dum fhidum] = min(abs(ffitIsf(idma:idmi) - thresh));
            TC.fhi{c}(p) = domIsf(fhidum+idma-1);  %high cutoff sfreq                      
            [dum flodum] = min(abs(ffitIsf(1:idma) - thresh));
            TC.flo{c}(p) = domIsf(flodum);  %low cutoff sfreq           
            
            TC.LPness{c}(p) = (ffitIsf(1)+phi(mi))/(ffitIsf(idma)+phi(mi));   
            
%             [ma idma] = max(ffitIsf); [mi idmi] = min(ffitIsf(idma:end));
%             thresh = ma*.5;
%             [dum fhidum] = min(abs(ffitIsf(idma:end) - thresh));
%             TC.fhi{c}(p) = domIsf(fhidum+idma-1);  %high cutoff sfreq
%                       
%             [dum flodum] = min(abs(ffitIsf(1:idma) - thresh));
%             TC.flo{c}(p) = domIsf(flodum);  %low cutoff sfreq
%             
%             TC.LPness{c}(p) = (ffitIsf(1))/(ffitIsf(idma));   
                       
            
            TC.Qfac{c}(p) = TC.sfpref{c}(p)/(TC.fhi{c}(p)-TC.flo{c}(p));
            TC.sfBW{c}(p) = log2(TC.fhi{c}(p)/TC.flo{c}(p)); %bandwidth in octaves    
            TC.sfBWLin{c}(p) = sfBWLin; %bandwidth in cyc/deg from sum of Gaussian fit  
            
            
            
            if TC.sfBW{c}(p) < .2  %high pass cells don't give values at zero, this is a reasonable cutoff                
                TC.sfBW{c}(p) = NaN;                
            end
            
            TC.LPness{c}(p) = (ffitIsf(1))/(ffitIsf(idma));           
            
            %TC.LPness{c}(p) = tcsf(1)/max(tcsf);   
            
            %Now get all the orientation metrics
            [TC.OMag{c}(p) TC.OAng{c}(p)] = orifind(tcoriI',DM.oridomI); 
            TC.opref{c}(p) = param(1);
            TC.orisig{c}(p) = param(2);
            TC.OAng{c}(p) = param(1);
            
            TC.RFsize{c}(p) = RFsize;
            
        else
            %spatial frequency:
            TC.sfpref{c}(p) = NaN;        
            TC.sfmag{c}(p) = 0;
            TC.Qfac{c}(p) = NaN;
            TC.sfBW{c}(p) = NaN;
            TC.sfBWLin{c}(p) = NaN;
            TC.LPness{c}(p) = NaN;
            TC.flo{c}(p) = NaN;
            TC.fhi{c}(p) = NaN;
            
            %orientation
            TC.OMag{c}(p) = NaN;     
            TC.OAng{c}(p) = NaN; 
            TC.opref{c}(p) = NaN;
            TC.orisig{c}(p) = NaN;
            
            TC.RFsize{c}(p) = NaN;
            
            TC.tcsfall_fit{c}(p,:) = NaN;
            TC.tcoriall_fit{c}(p,:) = NaN;
            
            
        end                
        
        %This is for the color tuning analyses;  It is redundant for each
        %cell 'c', but it is convenient later
        TC.tccolorall{1}(p,c) = sqrt(norm(TC.tcoriall{c}(p,:)) * norm(TC.tcsfall{c}(p,:)));
        
        
%         semilogx(DM.sfdom,tcsf,['.' colorid{c} '-']), hold on
%         plot(domIsf,ffitIsf,'k')
%         axis tight

                plot(DM.oridom,tcori,['.' colorid{c} '-']), hold on
                plot(DM.oridomI,ffit_ori{p},'k')
                axis tight

        %
        %         title(num2str(TC.Qfac{c}(p)))

%        imagesc(tcorisfraw)
%         hold on
%         contour(ffitorisf,(max(ffitorisf(:))+min(ffitorisf(:)))/2,'k')
         title(num2str(sigma))
         
         %Get stuff to plot the examples
         if ~isempty(find(p == idExamp))
            kern_examp{c}{idex} = tcorisfraw;   
            
            [idori idsf] = find(tcorisfraw == max(tcorisfraw(:)));
            tcoursema_examp{c}{idex} = squeeze(kernTplot(idori,idsf,:));            
            tsigcoursema_examp{c}{idex} = squeeze(kernsigTplot(idori,idsf,:));
            
            [idori idsf] = find(tcorisfraw == min(tcorisfraw(:)));
            tcoursemi_examp{c}{idex} = squeeze(kernTplot(idori,idsf,:));            
            tsigcoursemi_examp{c}{idex} = squeeze(kernsigTplot(idori,idsf,:));
            
            tcoriEx{c}{idex} = tcori;
            tcsfEx{c}{idex} = tcsf;
            
            ffitoriEx{c}{idex} = ffit_ori{p};
            ffitsfEx{c}{idex} = ffit_sf{p};
            
            idex = idex+1;
         end
        
    end

end
%This is dumb, but convenient later
TC.tccolorall{2} = TC.tccolorall{1};
TC.tccolorall{3} = TC.tccolorall{1};

if length(DM.colordom) > 1
    getTCfrom2colors
end

%% Compute F1/F0

sig = 50;
dom = DM.taudom-mean(DM.taudom);
psmooth = exp(-dom.^2/(2*sig^2));
psmooth = psmooth/sum(psmooth);
psmooth = ones(4,1)*psmooth;
psmooth = abs(fft(psmooth,[],2));

phasedomI = 0:359;

figure
for c = 1:length(DM.colordom)
    idex = 1;
    for p = 1:MK.Ncell
        subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
        

        [dum oriid] = min(abs(TC.OAng{1}(p)-DM.oridom));
        [dum sfid] = min(abs(TC.sfpref{1}(p)-DM.sfdom));
        
        kdum = kernC{c}{p}; domdum = DM.oridom;
        if length(DM.oridom) >= 50
            kdum = (kdum(1:2:end,:,:,:) + kdum(2:2:end,:,:,:))/2;
            domdum = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
        end
        [dum oriid] = min(abs(TC.OAng{c}(p)-domdum));

        kernplot = squeeze(kdum(oriid,sfid,:,:));  %phase and time
        kernplot = ifft(fft(kernplot,[],2).*psmooth,[],2);
        
        TC.tcphase{c}(p,:) = squeeze(kernplot(:,TC.tauID{c}(p)));  %use same time delay as ori/sf curves

        f1 = 2*mean(TC.tcphase{c}(p,:).*exp(1i*DM.phasedom*pi/180)); %peak amplitude
        f0 = mean(TC.tcphase{c}(p,:));

        TC.phase{c}(p) = angle(f1)*180/pi;       

        TC.F1F0{c}(p) = 2*abs(f1)/f0;  %peak-to-peak amplitude, over mean
        
        phasefit = abs(f1)*cos(phasedomI*pi/180 - TC.phase{c}(p)*pi/180) + f0;
        
        if TC.SNR{c}(p) < 1
            TC.phase{c}(p) = NaN;
            TC.F1F0{c}(p) = NaN;
        end

        %The Nishimoto et al version 
%         f1 = abs(TC.tcphase{c}(p,1)-TC.tcphase{c}(p,3)) + abs(TC.tcphase{c}(p,2)-TC.tcphase{c}(p,4));
%         f0 = mean(TC.tcphase{c}(p,:));
%         TC.F1F0{c}(p) = abs(f1)/f0;

        plot(DM.phasedom,TC.tcphase{c}(p,:),['.-' colorid{c}]), hold on
        plot(phasedomI,phasefit,'k'), axis tight

        title(num2str(TC.F1F0{c}(p)))

        %Get stuff to plot the examples
        if ~isempty(find(p == idExamp))

            tcphaseEx{c}{idex} = TC.tcphase{c}(p,:);
           
            ffitphaseEx{c}{idex} = phasefit;

            idex = idex+1;
        end

    end

end

figure,hist(TC.phase{1})

for c = 1:length(DM.colordom)
    id = find(isnan(TC.OMag{c}));
    TC.F1F0{c}(id) = NaN;
end


%id = find(TC.F1F0{c}<0);
%TC.F1F0{c}(id) = 0;
%id = find(TC.F1F0{c}>2);
%TC.F1F0{c}(id) = 2;
figure,hist(TC.F1F0{1},linspace(0,3,20))


%% Compute the spatial phase for a single ORI and SF
% clear oriid sfid
% for c = 1:length(DM.colordom)
%     popOpref = angle(nansum(exp(1i*TC.OAng{1}*2*pi/180)))/2;
%     if popOpref<0;
%         popOpref = popOpref + 180;
%     end
% 
%     id = find(~isnan(TC.sfpref{c}));
%     popSpref = geomean(TC.sfpref{c}(id));
% 
%     [dum oriid{c}] = min(abs(popOpref-DM.oridom));
%     [dum sfid{c}] = min(abs(popSpref-DM.sfdom));
% end
% 
% for p = 1:MK.Ncell
% 
%     for c = 1:length(DM.colordom)
%         
%         kernplot = squeeze(kernC{c}{p}(oriid{c},sfid{c},:,:));  %phase and time
%         kernplot = ifft(fft(kernplot,[],2).*psmooth,[],2);
%         
%         tcphase = squeeze(kernplot(:,TC.tauID{c}(p)));  %use same time delay as ori/sf curves
%         
%         tcphase = phi(tcphase);        
%         f1 = sum(tcphase.*exp(1i*DM.phasedom'*pi/180));
%         
%         TC.phase{c}(p) = angle(f1)*180/pi;
%         
%         %The Nishimoto et al version
% %         f1 = abs(tcphase(1)-tcphase(3)) + abs(tcphase(2)-tcphase(4));
% %         f0 = sum(tcphase);
% %         TC.F1F0{c}(p) = 2*f1/f0;
%         
%     end
% 
% end
% 
% for c = 1:length(DM.colordom)
%     id = find(isnan(TC.OMag{c}));
%     TC.phase{c}(id) = NaN;
% end

%% Compute 1D spatial profile

% sig = 50;
% dom = DM.taudom-mean(DM.taudom);
% psmooth = exp(-dom.^2/(2*sig^2));
% psmooth = ones(4,1)*psmooth;
% psmooth = abs(fft(psmooth,[],2));
% for i = 1:length(DM.sfdom)
%     smoother(i,:,:) = psmooth;
% end
% 
% figure
% for p = 1:MK.Ncell
%     subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
% 
%     for c = 1:length(DM.colordom)
%         
%         [dum oriid] = min(abs(TC.OAng{c}(p)-DM.oridom));
%         
%         kernplot = squeeze(kernC{c}{p}(oriid,:,:,:));  %spatial frequency, phase, and time
%         kernplot = ifft(fft(kernplot,[],3).*smoother,[],3);
%         
%         LineF = squeeze(kernplot(:,:,TC.tauID{c}(p)));  %use same time delay as ori/sf curves
%         
%         
%         
%         plot(DM.phasedom,tcphase,['.' colorid{c} '-']), hold on
%         
%         title(num2str(TC.F1F0{c}(p)))
%         
%     end
% 
% end

%% Plot examples
for i = 1:length(DM.sfdomI)
    sfdomcell{i} = num2str(num2str([round(DM.sfdomI(i)*10)/10]));
end

if ~isempty(idExamp)
    for c = 1:length(DM.colordom)
        figure
        Nex = length(kern_examp{1});
        for i = 1:Nex
            subplot(Nex*2,3,i*6-5)
            imagesc(kern_examp{c}{i})
            xl = str2num(get(gca,'XtickLabel'));
            set(gca,'XtickLabel',round(DM.sfdom(xl)*10)/10)
            yl = str2num(get(gca,'YtickLabel'));
            set(gca,'YtickLabel',round(DM.oridom(yl)*10)/10,'TickDir','out'), axis square

            subplot(Nex*2,3,i*6-4)

            fill([DM.taudom fliplr(DM.taudom)],[tcoursema_examp{c}{i}-tsigcoursema_examp{c}{i}; flipud(tcoursema_examp{c}{i}+tsigcoursema_examp{c}{i})]',[.0 .0 1])
            hold on
            plot(DM.taudom,tcoursema_examp{c}{i},'k'),xlabel('ms')
            %ylim([min(tcoursema_examp{c}{i}-tsigcoursema_examp{c}{i})-.1 max(tcoursema_examp{c}{i}+tsigcoursema_examp{c}{i})+.1])

            fill([DM.taudom fliplr(DM.taudom)],[tcoursemi_examp{c}{i}-tsigcoursemi_examp{c}{i}; flipud(tcoursemi_examp{c}{i}+tsigcoursemi_examp{c}{i})]',[1 .0 .0])
            hold on
            plot(DM.taudom,tcoursemi_examp{c}{i},'k'),xlabel('ms')
            xlim([DM.taudom(1) DM.taudom(end)])

            subplot(Nex*2,3,i*6-2)
            plot(DM.oridom,tcoriEx{c}{i},'.k'),xlabel('ori')
            hold on, plot(DM.oridomI,ffitoriEx{c}{i},'r')
            set(gca,'Xtick',[0 90 180]), xlim([-5 180])
            %ylim([min(ffitoriEx{c}{i})-.1 max(ffitoriEx{c}{i})+.1])

            subplot(Nex*2,3,i*6-1)
            plot(log2(DM.sfdom),tcsfEx{c}{i},'.k')
            hold on, plot(log2(DM.sfdomI),ffitsfEx{c}{i},'r'), %xlim([DM.sfdomI(1) DM.sfdomI(end)])
            %ylim([min(ffitsfEx{c}{i})-.1 max(ffitsfEx{c}{i})+.1])
            set(gca,'XTick',log2(DM.sfdom(1:2:end)))
            set(gca,'XTickLabel',round(DM.sfdom(1:2:end)*10)/10)
            xlim([log2(DM.sfdom(1)-.1) log(DM.sfdom(end))+1.5])

            subplot(Nex*2,3,i*6)
            plot([DM.phasedom 360],[tcphaseEx{c}{i} tcphaseEx{c}{i}(1)],'.k'),xlabel('phase')
            hold on, plot(phasedomI,ffitphaseEx{c}{i},'r')
            set(gca,'Xtick',[0 180 360]), xlim([-5 365])
            %ylim([min(ffitoriEx{c}{i})-.1 max(ffitoriEx{c}{i})+.1])

        end
    end
end


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
tcori = cir