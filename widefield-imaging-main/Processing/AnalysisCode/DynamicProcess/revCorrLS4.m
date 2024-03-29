function kern = revCorrLS4(cellMat,trialdom,hh)

%4 "deconvolves" the combined linear portion using a boot strap

%3 corrects for the nonlinearity w/ feedforward model, and then computes
%the ARMA model

global ACQinfo maskS Analyzer G_RChandles

DC = 1;
polyorder = 1;

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);
Ncell = length(nID);

ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
try
    load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')
catch
    load(['e:\2p_data\' Analyzer.M.anim '\log_files\' expt],'frate')
end

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 

[domains seqs] = getSeqInfo;

%%%%%%%%%%%%%%%%%%%

oridom = domains{1}.oridom;
sfdom = domains{1}.sfdom;

paramP = length(oridom)*length(sfdom);

tauP = round(1800/acqPeriod)+1;

%tauP = 40;
taudom = 0:acqPeriod:((tauP-1)*acqPeriod);
covMat = 0;
xCorr = cell(1,Ncell);
covMat = cell(1,Ncell);
for i = 1:Ncell
    xCorr{i} = 0;
    covMat{i} = 0;
end

figure
for p = 7:7

    pID = nID(p);
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    tauConst(p) = NaN;

    for trialid = 1:length(trialdom)

        T = trialdom(trialid);
        [cond rep] = getcondrep(T);
        
        hper = gethper(cond);         
        Tupdate = Tf*hper;   

        y = squeeze(cellMat{cond}(pID,:,rep));     

        %y = squeeze(mean(cellMat{cond}(:,:,rep),1));
        y = processTcourse(y,hh,polyorder,acqPeriod);
        
        tdom = (0:length(y)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain
        %of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;

        HHdum = zeros(length(y),length(oridom)*length(sfdom));
        pid = 1;
        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)

                id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf));
                %id = find(seqs{T}.oriseq == oridom(ori));

                if ~isempty(id)

                    stimes = (id-1)*Tupdate + getparam('predelay')*1000; %Stimulus times (ms)

                    for i = 1:length(stimes)
                        [dum idx1] = min(abs(tdom_pix - (stimes(i)+taudom(1)))); %find time sample that is closest to the beginning of response window
                        HHdum(idx1,pid) = 1;
                    end
                end
                pid = pid+1;
            end
        end
        %HHdum = zscore(HHdum);
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*tauP); %preallocate
        HHdum = [.01*randn(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdum(:,1))
            chunk = squeeze(HHdum(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end

        covMat{p} = covMat{p} + HH'*HH;

        xCorr{p} = HH'*y(:) + xCorr{p};

    end

    %covMat{p} = covMat{p}.*eye(size(covMat{p}));
    params{p} = (covMat{p})\xCorr{p};
    %params{p} = xCorr{p};
    
    paramsdum = params{p};

    if DC
        xCorr{p} = xCorr{p}(1:end-1);
        paramsdum = paramsdum(1:end-1);  %Get rid of DC shift
    end

    kerndum = reshape(paramsdum,paramP,tauP);

    for i = 1:length(kerndum(1,:))
        dum = kerndum(:,i);
        kern{p}(:,:,i) = reshape(dum,length(sfdom),length(oridom))';  %first dimension will be ori
    end
    
%This is the same thing as subtracting the response at time zero while
%building the kernel...
%     for i = 1:length(oridom)
%         for j = 1:length(sfdom)
%             kern{p}(i,j,:) = kern{p}(i,j,:) - kern{p}(i,j,1);
%         end
%     end
    

    %kernplot = fliplr(squeeze(mean(kern{p}(:,1:3,:),2)));  %Mean over sfreq... Get ori/time kernel
    kernplot = fliplr(squeeze(mean(kern{p}(3:6,:,:))));  %Get sf/time kernel
    %subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    subplot(2,2,1)
    imagesc(taudom,sfdom,kernplot)
%     set(gca,'YTickLabel',sfdom);
%     labs = str2num(get(gca,'XTickLabel'));
%     set(gca,'XTickLabel',round(labs*acqPeriod)/1000);
    title(num2str(tauConst(p)))
    
    [idy idx] = find(kernplot == max(kernplot(:)));
    subplot(2,2,2)
    plot(sfdom,kernplot(:,idx-1:idx+1)), legend('-T','peak','+T')
    subplot(2,2,3)
    plot(taudom,mean(kernplot(idy:idy,:),1))
     drawnow
    
    %imagesc(fliplr(Xcorrkern))

    %plot(kern(6,:))
end



%Now use bootstrap to fit the ARMA model


tauPARMA = round(200/acqPeriod)+1;

taudomARMA = 0:acqPeriod:((tauPARMA-1)*acqPeriod);
clear kern
covMat = 0;
xCorr = cell(1,Ncell);
covMat = cell(1,Ncell);
for i = 1:Ncell
    xCorr{i} = 0;
    covMat{i} = 0;
end

ARord = 1;

figure
for p = 7:7

    pID = nID(p);
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    tauConst(p) = NaN;

    for qq = 1:20         
        
        N = 5000;
        HHdum = rand(N,length(oridom)*length(sfdom));

        %Generate Bootstrap output
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*length(taudom)); 
        HHdumpad = [0*randn(tauP-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauP:length(HHdumpad(:,1))
            chunk = squeeze(HHdumpad(z-tauP+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauP+1,:) = chunk;
        end
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end
        
        y = zscore(HH*params{p});  %BS output 
        
        %Now fit the ARMA
        HH = zeros(length(HHdum(:,1)),length(HHdum(1,:))*length(taudomARMA)); %preallocate
        HHdumpad = [0*randn(tauPARMA-1,length(HHdum(1,:))); HHdum];  %Pad with zeros
        for z = tauPARMA:length(HHdumpad(:,1))
            chunk = squeeze(HHdumpad(z-tauPARMA+1:z,:))';
            chunk = chunk(:)';
            HH(z-tauPARMA+1,:) = chunk;
        end
        
        if DC
            HH = [HH ones(length(HH(:,1)),1)];
        end        
        
        for q = 1:ARord
            HH = [HH [zeros(q,1); y(1:end-q)]];
        end               

        covMat{p} = covMat{p} + HH'*HH;

        xCorr{p} = HH'*y(:) + xCorr{p};

    end

    %covMat{p} = covMat{p}.*eye(size(covMat{p}));
    params{p} = (covMat{p})\xCorr{p};
    %params{p} = xCorr{p};
    
    paramsdum = params{p};

    alpha = params{p}(end-ARord+1:end);
    tauConst(p) = acqPeriod/-log(alpha(1));
    paramsdum = params{p}(1:end-ARord);
    xCorr{p} = xCorr{p}(1:end-ARord);

    if DC
        xCorr{p} = xCorr{p}(1:end-1);
        paramsdum = paramsdum(1:end-1);  %Get rid of DC shift
    end

    kerndum = reshape(paramsdum,paramP,tauPARMA);
    
    for i = 1:length(kerndum(1,:))  %loop through each time point
        dum = kerndum(:,i);
        kern{p}(:,:,i) = reshape(dum,length(sfdom),length(oridom))';  %first dimension will be ori
    end
    
    %subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

    plotkerns(kern{p},oridom,sfdom,taudomARMA,acqPeriod,tauConst(p))
    
    %imagesc(fliplr(Xcorrkern))

    %plot(kern(6,:))
end

clear y
N = 80;
x = [zeros(1,ARord) 1 zeros(1,N-1)];
y = zeros(1,ARord);
for i = ARord+1:length(x)
    y(i) = x(i);
    for j = 1:ARord
        y(i) = y(i) + y(i-j)*alpha(j);
    end
end
tdom = (0:(length(y)-1))*acqPeriod;
tdom = tdom - acqPeriod*ARord;
figure,plot(tdom,y,'-ok')
[dum id] = min(abs(y-exp(-1)));
tau = tdom(id);
title(['tau = ' num2str(tau)])

function y = processTcourse(y,hh,polyorder,sp)

%mu = mean(y);
%First subtract a low-order polynomial fit:
yfit = polyfitter(y,polyorder);
y = y-yfit';

%Linear Bandpass Filter from GUI
if ~isempty(hh)
    y = ifft(fft(y).*hh);
end
%[y noise] = wiener2(y, [1 round(210/sp)])
%y = zscore(y);

%figure,plot(abs(fft(y)))

%Get rid of any peaks around the breathing rate
fdom = linspace(0,1/(sp/1000),length(y)+1);
fdom = fdom(1:end-1);
atten = [.2 .15 .1 .15 .2];
idbreath = find(fdom<2.5 & fdom>1);
Np = 0;
Npeaks = 3;
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
ydum = ifft(fft(y-mean(y)).*abs(fft(h)));
y = y - ordfilt2(ydum, 20, ones(1,600));

y = y-mean(y);

y = zscore(y);
%y = y.*sign(y);


function xfit = polyfitter(x,order)

dom = (0:length(x)-1)';

H = ones(length(dom),order+1);  %last column for DC
for i = 1:order   
    H(:,i) = dom.^i;    
end

p = inv(H'*H)*H'*x';
xfit = H*p;

function x = invertNL(x,y,yhat)


Nbins = 5;
ptsbin = floor(length(y_all{p})/Nbins);
[yhatS id] = sort(yhat_all{p});
yS = y_all{p}(id);
for i = 1:Nbins
    idsamp = ((i-1)*ptsbin+1):i*ptsbin;
    if i == Nbins
        idsamp = ((i-1)*ptsbin+1):length(y_all{p});
    end
    sampyhat = yhatS(idsamp);
    sampy = yS(idsamp);
    
    hatdum = sampyhat-trimmean(sampyhat,10);
    dum = sampy - trimmean(sampy,10);
    slpe = regress(dum,hatdum);
    
    
    muyhat(i) = trimmean(sampyhat,10);
    muy(i) = trimmean(sampy,10);
    sigyhat(i) = nanstd(sampyhat)/sqrt(length(sampyhat));
    sigy(i) = nanstd(sampy)/sqrt(length(sampy));
end


function plotkerns(kern,oridom,sfdom,taudom,acqPeriod,tauConst)

kernplot = fliplr(squeeze(mean(kern(:,1:3,:),2)));  %Mean over sfreq... Get ori/time kernel
figure
subplot(4,2,1)
imagesc(0:length(taudom)-1,0:length(oridom)-1,kernplot)
labs = str2num(get(gca,'YTickLabel'));
set(gca,'YTickLabel',labs*(oridom(2)-oridom(1)));
labs = str2num(get(gca,'XTickLabel'));
set(gca,'XTickLabel',round(labs*acqPeriod)/1000);
title(num2str(tauConst))

[idy idx] = find(kernplot == max(kernplot(:)));
subplot(4,2,2)
plot(oridom,kernplot(:,idx-1:idx+1)), xlim([oridom(1) oridom(end)]), %legend('-T','peak','+T')
subplot(4,2,3)
plot(taudom,mean(kernplot(idy:idy,:),1))

kernplot = fliplr(squeeze(mean(kern(3:6,:,:))));  %Get sf/time kernel
subplot(4,2,5)
imagesc(0:length(taudom)-1,0:length(sfdom)-1,kernplot)
labs = str2num(get(gca,'YTickLabel'));
set(gca,'YTickLabel',sfdom(1)+labs*(sfdom(2)-sfdom(1)));
labs = str2num(get(gca,'XTickLabel'));
set(gca,'XTickLabel',round(labs*acqPeriod)/1000);


[idy idx] = find(kernplot == max(kernplot(:)));
subplot(4,2,6)
plot(sfdom,kernplot(:,idx-1:idx+1)), xlim([sfdom(1) sfdom(end)]), %legend('-T','peak','+T')
subplot(4,2,7)
plot(taudom,mean(kernplot(idy:idy,:),1))
drawnow
    