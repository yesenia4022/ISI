function Grandposplots2(kern,tauN)


global ACQinfo G_RChandles cellS maskS DM MK Analyzer TC

%%%%

DM = struct; MK = struct; %Reset these

%Make/store another structure related to the maskS, 'MK'
MK.Ncell = length(cellS.kernAll);  

MK.masklabel = bwlabel(maskS.bwCell{1},4);
MK.celldom = unique(MK.masklabel);
[MK.nID] = getNeuronMask;

for p = 1:MK.Ncell
    [idcelly idcellx] = find(MK.masklabel == MK.celldom(MK.nID(p)));
    MK.CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

%Get/store the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string') ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with dtau spacing, and end at an estimate of kernDel(2)

%Get/store functional domains
trialdom = 1:1:getnotrials;
eval(['dT = ' get(G_RChandles.dropTrials,'string')])
trialdom(dT) = [];

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load([logfileroot Analyzer.M.anim '\' expt],'domains')

DM.oridom = domains.oridom;
DM.xdom = domains.xdom;
DM.colordom = domains.colordom;
DM.taudom = taudom;


%Get image dimensions
%xVolt = ACQinfo.scanAmplitudeX/ACQinfo.zoomFactor; %Scanimage 3.6
%yVolt = ACQinfo.scanAmplitudeY/ACQinfo.zoomFactor;
xVolt = 3/ACQinfo.zoomFactor;
yVolt = 3/ACQinfo.zoomFactor;

'WrongScanAmplitude'

xmic = 94*xVolt-2;  %Kristina fit these lines
ymic = 135*yVolt+0.5;
xmicperpix = xmic/ACQinfo.pixelsPerLine;
ymicperpix = ymic/ACQinfo.linesPerFrame;


 
%Downsample kernel and domains
if ~isempty(find(isnan(kern{1})))
    for p = 1:length(kern)
        
        dumA = kern{p}(:,1:2:end,:,:,:,:); dumB = kern{p}(:,2:2:end,:,:,:,:);
        id = find(isnan(dumA)); dumA(id) = dumB(id);
        id = find(isnan(dumB)); dumB(id) = dumA(id);

        kern{p} = (dumA + dumB)/2; %oridomain        

        dumA = kern{p}(1:2:end,:,:,:,:,:); dumB = kern{p}(2:2:end,:,:,:,:,:);
        id = find(isnan(dumA)); dumA(id) = dumB(id);
        id = find(isnan(dumB)); dumB(id) = dumA(id);

        kern{p} = (dumA + dumB)/2; %spatial domain
        
    end
    
    DM.oridom = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
    DM.xdom = (DM.xdom(1:2:end) + DM.xdom(2:2:end))/2;

    'Downsampling because not enough presentations'
end


for c = 1:length(DM.colordom)
    for p = 1:length(kern)
        kernC{c}{p} = squeeze(kern{p}(:,:,:,c,:));
    end    
end


ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  %it will start at exactly kernDel(1) with acqPeriod spacing, and end at an estimate of kernDel(2)

delayWin = [200 600];  %assume the peak response is within this time window
[dum delayWinID(1)] = min(abs(delayWin(1)-taudom));
[dum delayWinID(2)] = min(abs(delayWin(2)-taudom));

Ncell = length(kern);

%Get the receptive field
tN = length(taudom); oN = length(DM.oridom); xN = length(DM.xdom);
kernsmooth = getSmoother([.5 .7 1 .7 .5],[.5 1 .5],[.3 1 .3],tN,oN,xN);
orismoother = getSmoother([.1 .5 1 .5 .1],[.5 1 .5],[.4 1 .4],tN,oN,xN);
possmoother = getSmoother([.1 .3 .6 1 .6 .3 .1],[.2 1 .2],[.2 .5 1 .5 .2],tN,oN,xN);
%possmoother = getSmoother([.1 .5 .8 1 .8 .5 .1],[.1 .3 1 .3 .1],[.1 .3 1 .3 .1],tN,oN,xN);


param(5) = -36;

figure
for c = 1:length(DM.colordom)
    TC.xpos{c} = NaN*ones(1,Ncell);
    TC.ypos{c} = NaN*ones(1,Ncell);
    
    figure
    for p = 1:Ncell

        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

        kernplot_b = squeeze(kernC{c}{p}(:,:,1,:));
        kernplot_w = squeeze(kernC{c}{p}(:,:,2,:));

        kerndum_b = ifftn(fftn(kernplot_b).*abs(fftn(kernsmooth)));
        kerndum_b = kerndum_b(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kerndum_w = ifftn(fftn(kernplot_w).*abs(fftn(kernsmooth)));
        kerndum_w = kerndum_w(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window

        kernplot_b = smoothkern(kernplot_b,possmoother);
        kernplot_w = smoothkern(kernplot_w,possmoother);

        [ma idma] = max(kerndum_b+kerndum_w,[],3); %find maxima from smoothed version
        [bestoriid bestposid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestposid) + delayWinID(1) - 1;

        im_b = kernplot_b(:,:,tau);
        im_w = kernplot_w(:,:,tau);

        padder = zeros(length(im_b(:,1)),0);
        im_b = [padder im_b padder];
        im_w = [padder im_w padder];

        degperpix = getparam('x_size')/length(im_b(1,:));

        %im_b = phi(im_b-prctile(im_b(:),.1));
        %im_w = phi(im_w-prctile(im_w(:),.1));

        %     for i = 1:length(DM.oridom)
        %         tc = im_b(i,:) + im_b(i,:);
        %         [param ffit varaccount] = Gaussfit(1:length(tc),tc,0);
        %         imfit(:,i) = ffit;
        %     end

        imfit = im_w'+im_b';
        %imfit = phi(imfit-prctile(imfit(:),50));
        %imfit(:,end/2:end) = 0;

        RF = iradon(imfit,DM.oridom,'none','linear',1,length(imfit(:,1)));
        
        %RFw = iradon(im_w',DM.oridom,'none','linear',1,length(imfit(:,1)));
        %RFb = iradon(im_b',DM.oridom,'none','linear',1,length(imfit(:,1)));
        
        %RF = RFw/norm(RFw(:)) - RFb/norm(RFb(:)); 
        
        %RF = backproj(imfit,DM.oridom);
        %RF = phi(RF-prctile(RF(:),95));

        hh = zeros(size(RF));
        hh(1:5,1:5) = [.1 .3 1 .3 .1]'*[.1 .3 1 .3 .1];
        RF = ifft2(fft2(RF).*abs(fft2(hh)));

        %     [TC.xpos{c}(p) TC.ypos{c}(p)] = getmaxLoc(RF,30);  %returns 'degrees'

        otc = max(imfit);

        [dum ohat] = orifind(otc',DM.oridom);

        [param ffit varaccount] = Gaussfit2Drot(RF,ohat);
        TC.xpos{c}(p) = param(2)*degperpix;
        TC.ypos{c}(p) = param(1)*degperpix;
        xsize{c}(p) = param(3)*degperpix;
        ysize{c}(p) = param(4)*degperpix;

        if ysize{c}(p)>xsize{c}(p)
            param(5) = param(5)+90;
        else
            param(5) = param(5);
        end
        if param(5) < 0
            param(5) = 180+param(5);
        end
        if param(5) > 180
            param(5) = param(5)-180;
        end
        Opref{c}(p) = ohat;
        Opreffit{c}(p) = param(5);
        Omagfit{c}(p) = max([xsize{c}(p) ysize{c}(p)])/min([xsize{c}(p) ysize{c}(p)]);  %aspect ratio

        %     [yma xma] = find(RF == max(RF(:)));
        %     xdom = linspace(0,1,length(RF(1,:)))*getparam('x_size');
        %     ydom = linspace(0,1,length(RF(:,1)))*getparam('x_size');   %don't use 'y_size'
        %     xmarg = mean(RF);
        %     ymarg = mean(RF,2)';
        %     xmarg = RF(yma,:);
        %     ymarg = RF(:,xma)';
        %     xmarg = xmarg/sum(xmarg);
        %     ymarg = ymarg/sum(ymarg);
        %     TC.xpos{c}(p) = sum(xmarg.*xdom);
        %     TC.ypos{c}(p) = sum(ymarg.*ydom);


        xdom = (0:length(RF(1,:))-1)*degperpix; %degrees
        ydom = (0:length(RF(:,1))-1)*degperpix;   %don't use 'y_size'
        RF = RF/max(RF(:));
        if p == 1
            [yma xma] = find(RF == max(RF(:)));
            yran = yma-5:yma+5;
            xran = xma-5:xma+5;
            yran = 1:length(ydom);
            xran = 1:length(xdom);
        end
        dim = size(RF);
        RF = RF(yran,xran); ffit = ffit(yran,xran); ydom = ydom(yran); xdom = xdom(xran);
        
        TC.RF{c}{p} = RF;

        imagesc(xdom,ydom,RF), %title(num2str(param(5)))
        hold on
        contour(xdom,ydom,ffit,.6,'k')
        axis image

        [idy idx] = find(imfit' == max(imfit(:)));

        %plot([im_b(idy,:)' im_w(idy,:)'])
        
        %imagesc(im_b-im_w)
        
        %imagesc([imfit flipud(imfit)])
        drawnow
        axis off
        if p == 1
            title([num2str(ydom(end)-ydom(1)) 'degrees'])
        end

        %     contour(xdom,ydom,RF,.5,'k'), axis ij
        %     hold on
        % %     ylim([1.5 3.8])
        % %     xlim([.2 2.3])
        %     drawnow

        %plot([mean(kernplot_b(4:5,:))' mean(kernplot_w(4:5,:))'])

    end
end

% thet = 40;
% xposp = TC.xpos*cos(thet) + TC.ypos*sin(thet);
% yposp = -TC.xpos*sin(thet) + TC.ypos*cos(thet);


    
figure
for c = 1:length(DM.colordom)
    TC.tcoriall{c} = [];
    TC.OAng{c} = [];
    for p = 1:length(kern)
    
        kernplot_b = squeeze(kernC{c}{p}(:,:,1,:));
        kernplot_w = squeeze(kernC{c}{p}(:,:,2,:));

        kerndum_b = ifftn(fftn(kernplot_b).*abs(fftn(kernsmooth)));
        kerndum_b = kerndum_b(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window
        kerndum_w = ifftn(fftn(kernplot_w).*abs(fftn(kernsmooth)));
        kerndum_w = kerndum_w(:,:,delayWinID(1):delayWinID(2));  %limit it to a reasonable time window

        kernplot_b = ifftn(fftn(kernplot_b).*abs(fftn(possmoother)));
        kernplot_w = ifftn(fftn(kernplot_w).*abs(fftn(possmoother)));

        [ma idma] = max(kerndum_b+kerndum_w,[],3); %find maxima from smoothed version
        [bestoriid bestposid] = find(ma == max(ma(:)));
        tau = idma(bestoriid,bestposid) + delayWinID(1) - 1;

        %     tcori_b = kernplot_b(:,bestposid,tau);  %Can't do this dumb assssss!!!
        %     tcori_w = kernplot_w(:,bestposid,tau);

        %     tcori_b = mean(kernplot_b(:,:,tau),2);  %Shouldn't do this either. If it is actually a projection then they will all be equal
        %     tcori_w = mean(kernplot_w(:,:,tau),2);
        %     tcori = tcori_b+tcori_w;

        tcori = max(kernplot_b(:,:,tau)+kernplot_w(:,:,tau),[],2);
        
        TC.tcoriall{c}(p,:) = tcori;

        [TC.OMag{c}(p) TC.OAng{c}(p)] = orifind(tcori,DM.oridom);

        [param ffit varacc] = Gaussfit(DM.oridom,tcori',1);

        if varacc < .6
            param = param*NaN;
        end

        TC.OAng{c}(p) = param(1);

        subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)

        plot(tcori)

    end

end



function smoother = getSmoother(ktau,kori,kpos,tN,oN,pN)

if oN == 1
    kori = 1;
end

ktau = [ktau zeros(1,tN-length(ktau))];
kpos = [kpos zeros(1,pN-length(kpos))];
kori = [kori zeros(1,oN-length(kori))];
kdum = kori'*kpos;
for i = 1:length(ktau)
    smoother(:,:,i) = kdum*ktau(i);
end
smoother = smoother/sum(smoother(:));

function kern = smoothkern(kern,smoother)

oriH = [.5 1 .5];
posH = [.5 1 .5]; Po = floor(length(posH)/2);
tauH = [.2 .5 1 .5 .2];

dim = [size(kern,1)*2 size(kern,2)+2*Po size(kern,3)];

kpad = zeros(dim);
for i = 1:size(kern,3)
    kdum = kern(:,:,i);
    kdum = [kdum; fliplr(kdum)];
    kdum = [kdum(:,1)*ones(1,Po) kdum kdum(:,end)*ones(1,Po)];
    kpad(:,:,i) = kdum;    
end

possmoother = getSmoother(tauH,oriH,posH,dim(3),dim(1),dim(2));
kern = ifftn(fftn(kpad).*abs(fftn(possmoother)));

kern = kern(1:dim(1)/2,Po+1:dim(2)-Po,:);

function [OMag OAng] = orifind(G,oridomain)

R = sum(G'.*exp(1i*oridomain*pi/90));
OAng = angle(R);                 %-pi to pi
OAng = OAng + pi*(1-sign(OAng+eps));  %0 to 2pi
OAng = OAng*90/pi;               %0 to 180
OMag = abs(R)/(sum(G));


function [xpos ypos] = getmaxLoc(RF,D)

xdom = linspace(0,1,length(RF(1,:)))*getparam('x_size'); %degrees
ydom = linspace(0,1,length(RF(:,1)))*getparam('x_size');   %don't use 'y_size'
[xmat ymat] = meshgrid(xdom,ydom);

xdomI = linspace(xdom(1),xdom(end),D*length(xdom));
ydomI = linspace(ydom(1),ydom(end),D*length(ydom));
[xmatI ymatI] = meshgrid(xdomI,ydomI);
RFI = interp2(xmat,ymat,RF,xmatI,ymatI);

[yma xma] = find(RFI == max(RFI(:)));
xpos = xdomI(xma);
ypos = ydomI(yma);


function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;