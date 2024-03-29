function kern = revCorrLS(cellMat,trialdom,hh)

global ACQinfo maskS Analyzer G_RChandles

ARflag = 1;
polyorder = 1;

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);
Ncell = length(nID);

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)

hper = getparam('h_per');

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 
Tupdate = Tf*hper;

[domains seqs] = getSeqInfo;

%%%%%%%%%%%%%%%%%%%

oridom = domains{1}.oridom;
sfdom = domains{1}.sfdom;

paramP = length(oridom);
tauP = 10;
covMat = 0;
xCorr = cell(1,Ncell);
covMat = cell(1,Ncell);
for i = 1:Ncell
    xCorr{i} = 0;
    covMat{i} = 0;
end

figure
for p = 1:7

    pID = nID(p);
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    for trialid = 1:length(trialdom)

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);

        y = squeeze(cellMat{cond}(pID,:,rep));        
        y = processTcourse(y,hh,polyorder,acqPeriod);
        
        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;

        HHdum = zeros(length(y),length(oridom));
        for ori = 1:length(oridom)
            %for sf = 1:length(sfdom)

            %id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
            id = find(seqs{T}.oriseq == oridom(ori));

            if ~isempty(id)

                stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                for i = 1:length(stimes)

                    [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                    HHdum(idx1,ori) = 1;
                end
            end
            %end
        end
        %HHdum = zscore(HHdum);
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [zeros(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        
        %HH = [HH ones(length(HH(:,1)),1)];
        if ARflag
            HH = [HH [0; y(1:end-1)']];
        end

        covMat{p} = covMat{p} + HH'*HH;

        xCorr{p} = HH'*y(:) + xCorr{p};

    end


    params{p} = inv(covMat{p})*xCorr{p};
    tauConst{p} = 0;
    paramsdum = params{p};
    if ARflag
        alpha = params{p}(end);
        tauConst{p} = acqPeriod/-log(alpha);
        paramsdum = params{p}(1:end-1);
        xCorr{p} = xCorr{p}(1:end-1);
    end
    %xCorr{p} = xCorr{p}(1:end-1);
    %paramsdum = paramsdum(1:end-1);  %Get rid of DC shift

    kern{p} = reshape(paramsdum,paramP,tauP);

    Xcorrkern = reshape(xCorr{p},paramP,tauP);

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    imagesc(fliplr(kern{p}))
    title(num2str(tauConst{p}))
    drawnow
    %imagesc(fliplr(Xcorrkern))

    %plot(kern(6,:))
end



%%Now Prediction



for p = 7:7

    pID = nID(p);
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    for trialid = 11:11
        T = trialdom(trialid);
        [cond rep] = getcondrep(T);

        y = squeeze(cellMat{cond}(pID,:,rep));
        y = processTcourse(y,hh,polyorder,acqPeriod);

        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;

        HHdum = zeros(length(y),length(oridom));
        for ori = 1:length(oridom)
            %for sf = 1:length(sfdom)

            %id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
            id = find(seqs{T}.oriseq == oridom(ori));

            if ~isempty(id)

                stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                for i = 1:length(stimes)

                    [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                    HHdum(idx1:idx1,ori) = 1;
                end
            end
            %end
        end

        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [zeros(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        
        %HH = [HH ones(length(HH(:,1)),1)];
        if ARflag
            MA = HH*params{p}(1:end-1);  %The moving average
            yhat = 0;
            for i = 1:length(MA);  %Now the autoregression
                yhat(i+1) = yhat(i)*params{p}(end) + MA(i);
            end
            yhat = yhat(2:end);
        else
            yhat = HH*params{p};
        end

        yhat = (yhat(50:end-100));
        y = (y(50:end-100));
        tdom = (0:length(y)-1)*acqPeriod/1000;
        
        figure
        plot(y)
        hold on
        plot(yhat,'r')
        title(num2str(tauConst{p}))
        drawnow
        
        Nbins = 10;
        ptsbin = floor(length(y(:))/Nbins);
        [yhatS id] = sort(yhat(:));
        yS = y(id);
        for i = 1:Nbins
            idsamp = ((i-1)*ptsbin+1):i*ptsbin;
            if i == Nbins
                idsamp = ((i-1)*ptsbin+1):length(y(:));
            end
            samplo = yhatS(idsamp);
            samphi = yS(idsamp);

            mulo(i) = trimmean(samplo,10);
            muhi(i) = trimmean(samphi,10);
            siglo(i) = nanstd(samplo)/sqrt(length(samplo));
            sighi(i) = nanstd(samphi)/sqrt(length(samplo));
        end
        figure,
        scatter(yhat,y,'.')
        
        hold on
        errorbar(mulo,muhi,sighi,'k','linewidth',2)
        
        corrcoef(y(:),yhat(:))
    end
end

function y = processTcourse(y,hh,polyorder,sp)

%mu = mean(y);
%First subtract a low-order polynomial fit:
yfit = polyfitter(y,polyorder);
y = y-yfit';

%Linear Bandpass Filter from GUI
if ~isempty(hh)
    y = ifft(fft(y).*hh);
end
%y = zscore(y);

%figure,plot(abs(fft(y)))

%Get rid of any peaks around the breathing rate
fdom = linspace(0,1/(sp/1000),length(y)+1);
fdom = fdom(1:end-1);
atten = [.2 .15 .1 .15 .2];
idbreath = find(fdom<2.5 & fdom>1);
Np = 0;
Npeaks = 0;
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
%y = y.*sign(y);


function xfit = polyfitter(x,order)

dom = (0:length(x)-1)';

H = [];
for i = 1:order   
    H = [H dom.^i];    
end
H = [H ones(length(H(:,1)),1)];

p = inv(H'*H)*H'*x';
xfit = H*p;