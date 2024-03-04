function make1DRF(kern,bwCell1,taudom)

global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2});

dim = [length(domains.oridom) length(domains.sfdom) length(domains.phasedom) length(domains.colordom)];

for i = 1:length(rseeds)
    
    eval(['colorS = rseed' num2str(i) '.colorseq;']);
    eval(['phaseS = rseed' num2str(i) '.phaseseq;']);
    eval(['sfS = rseed' num2str(i) '.sfseq;']);
    eval(['oriS = rseed' num2str(i) '.oriseq;']);
    
    colorseq{i} = domains.colordom(colorS);   
    phaseseq{i} =  domains.phasedom(phaseS);  
    sfseq{i} =  domains.sfdom(sfS);     
    oriseq{i} =  domains.oridom(oriS);  
 
end

% if round(oridom(end)) == 999
%     oridom = oridom(1:end-1);
% end

%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;

Ncell = length(kern);
NT = getnotrials;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make smoothing functions (used to find maxima)
ktau = [.1 .5 1 .5 .1];
ktau = [ktau zeros(1,length(taudom)-length(ktau))];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot sf curves

stimsize = getparam('x_size');
degdom = linspace(0,stimsize,100);

figure
for p = 1:Ncell    

    kernplot = kern{p};

    kernplot = ifftn(fftn(kernplot).*abs(fftn(kernsmooth)));
    
    kernplot = kernplot-mean(kernplot(:));


    %kernplot = kern{p};    
    Frespdum = squeeze(mean(kernplot(:,:,:,2:4),4));
    Fresp(:,:,1) = Frespdum(:,:,1) - Frespdum(:,:,3); %combine 0 and 180
    Fresp(:,:,2) = Frespdum(:,:,2) - Frespdum(:,:,4); %combine 90 and 270
    
    Freal = Fresp(:,:,1);
    Fimag = Fresp(:,:,2);
    F = Freal - 1i*Fimag;

    for ori = 1:length(oridom)               
        
        RF = 0;
        for sf = 1:length(sfdom)               
                %RF = Freal(ori,sf) * cos(degdom*sfdom(sf)*2*pi) + RF;
                %RF = Fimag(ori,sf) * sin(degdom*sfdom(sf)*2*pi) + RF;
                RF = F(ori,sf)*exp(1i*(degdom*sfdom(sf)*2*pi)) + RF;
        end
        for sf = 1:length(sfdom)             
                RF = F(ori,sf)*exp(1i*(degdom*-sfdom(sf)*2*pi)) + RF;
        end
        

        RFradon(:,ori) = RF(:);
        
    end

    RF = iradon((RFradon),oridom);
    
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    imagesc(abs(RF))
    

end

