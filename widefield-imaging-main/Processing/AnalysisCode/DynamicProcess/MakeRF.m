function RF = MakeRF(kern,Ntau)

global Analyzer


%%%%

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])

%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;
colordom = domains.colordom;

Ncell = length(kern);
NT = getnotrials;

Ntau = Ntau+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)
ktau = [.1 .5 1 .5 .1];
ktau = [ktau zeros(1,Ntau-length(ktau))];
ksf = [.1 1 .1];
ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
kori = [.5 1 .5];
kori = [kori zeros(1,length(oridom)-length(kori))];
kphase = [.05 1 .05];
kphase = [kphase zeros(1,length(phasedom)-length(kphase))];

kdum = kori'*ksf;
for i = 1:length(ktau)
    kerndum(:,:,i) = kdum*ktau(i);
end
for i = 1:length(kphase)
    kernsmooth(:,:,i,:) = kerndum*kphase(i);
end
kernsmooth = kernsmooth/sum(kernsmooth(:));


for i = 1:length(colordom)
    for p = 1:length(kern)
        kernC{i}{p} = squeeze(kern{p}(:,:,:,i,:));
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot sf curves

stimsize = getparam('x_size');

Wcm = stimsize/360*(2*pi*Analyzer.M.screenDist);
pixpercm = Analyzer.M.xpixels/Analyzer.M.screenXcm;
Npix = Wcm*pixpercm;

stimsize = 1/(sfdom(1))*2/3;
phi = linspace(-stimsize/2,stimsize/2,Npix);

[x y] = meshgrid(phi,phi);

figure
for p = 1:Ncell
    
    kernplot = kernC{1}{p};
    
    %kernplot = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));    
    
    Frespdum = squeeze(mean(kernplot(:,:,:,3:5),4));
    %Frespdum = Frespdum - min(Frespdum(:));
    Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2)+Frespdum(:,:,3)+Frespdum(:,:,4));
    %Frespdum = squeeze(Frespdum(:,:,1)-Frespdum(:,:,3));
    
    RF = 0;
    for ori = 1:length(oridom)
    %for ori = 1:8
        %for sf = 1:length(sfdom)
        for sf = 1:5
            
            xphi = sfdom(sf)*x*2*pi; %radians
            yphi = sfdom(sf)*y*2*pi;
            
            xp = xphi*cos(oridom(ori)*pi/180) + yphi*sin(oridom(ori)*pi/180); %radians
            grat = Frespdum(ori,sf)*cos(xp - pi/4);  %Hartley inverse trx
            RF = RF + grat;
            
            %xp = xphi(1,:);
            %RFiradom(:,ori) = Frespdum(ori,sf)*cos(xp-pi/4);
            
            
        end
    end

    %tc = mean(Frespdum,3);
    tc = squeeze(mean(Frespdum(:,2:4),2));
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    %imagesc(RFiradom)

    imagesc(RF)
    axis off
    
end