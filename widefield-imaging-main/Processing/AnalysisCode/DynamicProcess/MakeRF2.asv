function RF = MakeRF2(kern,Ntau)

global Analyzer DM cellS

%%%%%%%%%%%%%%%%%%%
oridom = DM.oridom;
sfdom = DM.sfdom;
phasedom = DM.phasedom;
colordom = DM.colordom;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)
% ktau = [.1 .5 1 .5 .1];
% ktau = [ktau zeros(1,Ntau-length(ktau))];
% ksf = [.1 1 .1];
% ksf = [ksf zeros(1,length(sfdom)-length(ksf))];
% kori = [.5 1 .5];
% kori = [kori zeros(1,length(oridom)-length(kori))];
% kphase = [.05 1 .05];
% kphase = [kphase zeros(1,length(phasedom)-length(kphase))];
% 
% kdum = kori'*ksf;
% for i = 1:length(ktau)
%     kerndum(:,:,i) = kdum*ktau(i);
% end
% for i = 1:length(kphase)
%     kernsmooth(:,:,i,:) = kerndum*kphase(i);
% end
% kernsmooth = kernsmooth/sum(kernsmooth(:));


for i = 1:length(colordom)
    for p = 1:length(kern)
        kernC{i}{p} = squeeze(kern{p}(:,:,:,i,:));
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ncell = length(cellS.muTime);

%plot sf curves

stimsize = getparam('x_size');

Wcm = stimsize/360*(2*pi*Analyzer.M.screenDist);
pixpercm = Analyzer.M.xpixels/Analyzer.M.screenXcm;
Npix = Wcm*pixpercm;

%stimsize = 1/(sfdom(1))*2/3;
phi = linspace(-stimsize/2,stimsize/2,Npix);

colors = {'r','g','b'};

Dec = 1;
RFall = 0
figure
for p = 1:Ncell
    
    %for c = 1:3
        kernplot = kernC{1}{p};
        
        %kernplot = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
        
        Frespdum = squeeze(mean(kernplot(:,:,:,9:15),4));
        %Frespdum = Frespdum - min(Frespdum(:));
        %Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2)-Frespdum(:,:,3)-Frespdum(:,:,4));
        %Frespdum = squeeze(Frespdum(:,:,1)+Frespdum(:,:,2));
        
        for ori = 1:length(oridom)
            
            RF = 0;
            for sf = 1:3
                
                for phase = 1:2
                    
                    xphi = sfdom(sf)*phi*2*pi/Dec; %radians
                    
                    RF = RF + Frespdum(ori,sf,phase)*sin(xphi-phasedom(phase)*pi/180);
                    
                end
                
                for phase = 3:4
                    
                    xphi = sfdom(sf)*phi*2*pi/Dec; %radians
                    
                    RF = RF - Frespdum(ori,sf,phase)*sin(xphi-phasedom(phase)*pi/180);
                    
                end
                
            end
            RFiradon(:,ori) = RF;
        end
        
        %tc = mean(Frespdum,3);
        tc = squeeze(mean(Frespdum(:,2:4),2));
        subplot(ceil(sqrt(Ncell/2)),ceil(sqrt(Ncell/2)),p)
        
        [dum id] = max(var(RFiradon));
        plot(RFiradon(:,id))
        ylim([-2 2])
        
        RF = iradon(RFiradon,oridom);
        %hold on
        %imagesc(RFiradon)
        
        %plot(mean(RFiradon,2),colors{c})
        
        %plot(squeeze(Frespdum(2,2,:)))
        
        RFall = RFall+RF;
        
        
        imagesc(RFiradon)
        %axis off
    %end
    
end

%figure,imagesc(RFall)