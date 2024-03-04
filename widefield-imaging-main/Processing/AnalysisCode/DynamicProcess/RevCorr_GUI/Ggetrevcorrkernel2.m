function Ggetrevcorrkernel2(cellMat,trialdom,hh)

%Ian Nauhaus

%Keep cellMat as an input because it filters it, and I don't want to
%create a new variable of the same size

%2 computes the standard deviation as well. 

global ACQinfo Analyzer cellS maskS G_RChandles

blankNorm = get(G_RChandles.blankNorm,'value');


%%%Get rid of the glia:

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)
[dum id0] = min(abs(taudom-0));  
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt],'frate')

if ~exist('frate')
    frate = 60.0044;
end

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file)

[domains seqs] = getSeqInfo(trialdom);


%%%%%%%%%%%%%%%%%%%

oridom = domains{trialdom(1)}.oridom;
sfdom = domains{trialdom(1)}.sfdom;
phasedom = domains{trialdom(1)}.phasedom;
colordom = domains{trialdom(1)}.colordom;

Ncell = length(nID);

cellS.muTime = cell(1,Ncell);
cellS.sigTime = cell(1,Ncell);
cellS.muBase = zeros(1,Ncell);
cellS.sigBase = zeros(1,Ncell);
cellS.muTimeBlank = cell(1,Ncell);
cellS.sigTimeBlank = cell(1,Ncell);


syncDiodeFlag = 1;

%Preallocate
Precountmat = getGratingPresentationCounts(trialdom);
for ori = 1:length(oridom)
    for sf = 1:length(sfdom)
        for phase = 1:length(phasedom)
            for color = 1:length(colordom)
                kernPre{ori,sf,phase,color} = zeros(Precountmat(ori,sf,phase,color),length(taudom));
            end
        end
    end
end

NT = getnotrials;
tcoursedum = 0;
figure
for p = 1:Ncell
tic
    pID = nID(p);

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom)); 
    countmatblank{p} = 0;
    
    %kern{p} = NaN*zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));    
    %kernblank{p} = zeros(1,length(taudom));


    for trialid = 1:length(trialdom)

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);

        hper = gethper(cond);
        Tupdate = Tf*hper;

        %         dsync = diff(cellS.synctimes{cond,rep});
        %         if any(abs(dsync-Tupdate/1000)>.100)
        %             'Warning: syncs may be messed up'
        %         end

        tcourse = squeeze(cellMat{cond}(pID,:,rep));
        
        id = find(isnan(tcourse));
        tcourse(id) = nanmean(tcourse);

        if ~isnan(sum(tcourse))  %some cells go out of the field of view

            tcourse = processTcourse(tcourse,hh,1,acqPeriod);

            %noise = fftshift(xcov(tcourse,'unbiased'));
            %noise = noise(end) - noise(end-1); %estimate of noise variance; assuming its white, additive, and independent from everything else

            %[tcourse noise] = wiener2(tcourse, [1 round(700/acqPeriod)],noise);

            %tcourse = myBigWeiner(1/500,acqPeriod,tcourse);

            %tcourse = myBigWeiner(1/500,acqPeriod,tcourse);

            tdom = (0:length(tcourse)-1)*acqPeriod;
            %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
            tdom_pix = tdom + tau_xy;

            %Deconvolution:
            %Ctrx = exp(-tdom/900);
            %tcourse = ifft(fft(tcourse)./fft(Ctrx));


            %idbase = find(tdom_pix<cellS.synctimes{cond,rep}(1)*1000 | tdom_pix>(cellS.synctimes{cond,rep}(end)*1000+500));
            %bLine = mean(tcourse(idbase));

            %tcourse = (tcourse-bLine)/bLine;

            if syncDiodeFlag
                dstimes = diff(cellS.synctimes{cond,rep}*1000);
                if find(dstimes>Tupdate*1.5 | dstimes<Tupdate*.5)
                    'messed up synctimes!!!!!!!!!'
                end
            end

            
            for ori = 1:length(oridom)
                for sf = 1:length(sfdom)
                    for phase = 1:length(phasedom)
                        for color = 1:length(colordom)

                            id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf) & seqs{T}.phaseseq == phasedom(phase)& seqs{T}.colorseq == colordom(color));

                            if ~isempty(id)

                                if syncDiodeFlag
                                    stimes = cellS.synctimes{cond,rep}(id)*1000;
                                else
                                    stimesdum = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)
                                end
                                
                                %stimesW = [30000 60000];
                                %stimes(find(stimes>stimesW(2) | stimes<stimesW(1))) = [];

                                for i = 1:length(stimes)

                                    [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                                    %dphs = tdom_pix(idx1) - (stimes(i)+taudom(1));  %Phase shift of the time window

                                    if ~isempty(idx1)
                                        idx1 = idx1(1);
                                        tpiece = idx1:idx1+length(taudom)-1;
                                        tbase = (idx1+id0-1)-1:(idx1+id0-1)-1;  %Use time zero as base

                                        if tpiece(1)>0 && tpiece(end)<length(tcourse)
                                            base = mean(squeeze(tcourse(tbase)));
                                            tcoursePiece = squeeze(tcourse(tpiece))';

                                            %tcoursePiece = (tcoursePiece-base);
                                            kern{p}{ori,sf,phase,color}(countmat{p}(ori,sf,phase,color)+1,:) = tcoursePiece;
                                            countmat{p}(ori,sf,phase,color) = countmat{p}(ori,sf,phase,color) + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end



            id = find(isnan(seqs{T}.oriseq));
            stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times
            %stimes = cellS.synctimes{cond,rep}(id)*1000;

            %stimes(find(stimes>stimesW(2) | stimes<stimesW(1))) = [];


            for i = 1:length(stimes)

                [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                if ~isempty(idx1)

                    idx1 = idx1(1);
                    tpiece = idx1:idx1+length(taudom)-1;
                    tbase = (idx1+id0-1)-1:(idx1+id0-1)-1;  %Use time zero as base

                    if tpiece(1)>0 & tpiece(end)<length(tcourse)
                        base = mean(squeeze(tcourse(tbase)));
                        tcoursePiece = squeeze(tcourse(tpiece))';

                        %tcoursePiece = (tcoursePiece-base);
                        kernblank{p}(countmatblank{p}+1,:) = tcoursePiece;
                        countmatblank{p} = countmatblank{p} + 1;
                    end
                end
            end

        end
    end

    id = find(countmat{p} == 0);
    if ~isempty(id)
        'Some stimuli never shown.'
    end

    if sum(countmatblank{1})
        cellS.muTimeBlank{p} = mean(kernblank{p});
        cellS.sigTimeBlank{p} = std(kernblank{p})/sqrt(length(kernblank{p})); %standard error
    end

    cellS.kernAll = kern;

    kern = [];

    for ori = 1:length(oridom)
        for sf = 1:length(sfdom)
            for phase = 1:length(phasedom)
                for color = 1:length(colordom)
                    if ~isempty(cellS.kernAll{p}{ori,sf,phase,color}) %if it didn't show this grating
                        N = length(cellS.kernAll{p}{ori,sf,phase,color}(:,1)); %Number of presentations

                        if blankNorm
                            cellS.muTime{p}(ori,sf,phase,color,:) = squeeze(mean(cellS.kernAll{p}{ori,sf,phase,color},1) - cellS.muTimeBlank{p});% / mean(cellS.muTimeBlank{p});
                            cellS.sigTime{p}(ori,sf,phase,color,:) = std(cellS.kernAll{p}{ori,sf,phase,color},[],1)/sqrt(N);% / mean(cellS.muTimeBlank{p});
                        else
                            cellS.muTime{p}(ori,sf,phase,color,:) = squeeze(mean(cellS.kernAll{p}{ori,sf,phase,color},1));
                            cellS.sigTime{p}(ori,sf,phase,color,:) = std(cellS.kernAll{p}{ori,sf,phase,color},[],1)/sqrt(N);
                        end
                    end
                end
            end
        end
    end

    %Subtract mean from each time point
    % for z = 1:length(cellS.muTime{p}(1,1,1,1,:))
    %
    %     dum = cellS.muTime{p}(:,:,:,:,z);
    %     cellS.muTime{p}(:,:,:,:,z) = cellS.muTime{p}(:,:,:,:,z) - mean(dum(:));
    %
    % end

    %Base is used for data selection later...
    %     dum = cellS.muTime{p}(:,:,:,:,1);
    %     cellS.muBase(p) = mean(dum(:));  %average response at time zero
    %     dum = cellS.sigTime{p}(:,:,:,:,1);
    %     cellS.sigBase(p) = mean(dum(:));  %average sigma response at time zero


    kernplot = cellS.muTime{p}(:,:,:,:,:);
    kernplot = mean(kernplot,3);  %mean across phase
    kernplot = mean(kernplot,4); %mean across color

    %kernblank = kern{p}(end,:);

    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    if size(kernplot,1) == 1 %if there was only one orientation
        kernplot = squeeze(nanmean(kernplot,1));  %mean across orientation
    else
        kernplot = squeeze(nanmean(kernplot,2));  %mean across spatial freqency
    end
    
    %blankmat = ones(length(kernplot(:,1)),1)*kernblank{p}(:)';

    %kernplot = (kernplot-blankmat);
    %kernplot = [kernplot; kernblank{p}(:)'];

    imagesc(round(taudom),sfdom,real(kernplot))
    %plot(cellS.muTimeBlank{p})

    drawnow
    toc
end


function y = processTcourse(y,hh,polyorder,sp)

id = find(isnan(y));
y(id) = nanmean(y);

%mu = mean(y);
%First subtract a low-order polynomial fit:
yfit = polyfitter(y,polyorder);
y = y-yfit';

%Linear Bandpass Filter from GUI
if ~isempty(hh)
    y = ifft(fft(y).*hh);
end
%[y noise] = wiener2(y, [1 round(300/sp)]);
%y = zscore(y);

%figure,plot(abs(fft(y)))

%Get rid of any peaks around the breathing rate
fdom = linspace(0,1/(sp/1000),length(y)+1);
fdom = fdom(1:end-1);
atten = [.4 .15 .1 .15 .4];
idbreath = find(fdom<2.5 & fdom>.5);
Np = 0;
Npeaks = 2;
for i = 1:Npeaks

    yf = fft(y);
    [ma id] = max(abs(yf(idbreath)));
    if ma > 3*std(abs(yf(idbreath))) + mean(abs(yf(idbreath)));
        idpeak = id+idbreath(1)-1;
        yf(idpeak-2:idpeak+2) = yf(idpeak-2:idpeak+2).*atten;

        yf = fliplr(yf);
        [dum id] = max(abs(yf(idbreath)));
        idpeak = id+idbreath(1)-1;
        yf(idpeak-2:idpeak+2) = yf(idpeak-2:idpeak+2).*atten;
        yf = fliplr(yf);
        y = real(ifft(yf));
        Np = Np+1;
    else
        break
    end

end

%hold on,
%plot(abs(fft(y)),'r')
%title(num2str(Np))

%Nonlinear highpass filter to subtract estimated baseline
h = fspecial('gaussian', [1 length(y)], 10); %Use stats from heavily smoothed version
ydum = ifft(fft(y).*abs(fft(h)));
y = y - ordfilt2(ydum, 20, ones(1,600));

id = find(y<prctile(y,50) & y>prctile(y,0));
y = y-median(y); 
y = y/std(y(id));

%y = zscore(y);
%y = y.*sign(y);


function xfit = polyfitter(x,order)

dom = (0:length(x)-1)';

H = ones(length(dom),order+1);  %last column for DC
for i = 1:order
    H(:,i) = dom.^i;
end

p = inv(H'*H)*H'*x';
xfit = H*p;