function XcorrPop

global maskS cellS ACQinfo

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 

nID = getNeuronMask;  %get the index values for the neurons
masklabel = bwlabel(maskS.neuronmask,4);
celldom = unique(masklabel);
celldom = celldom(2:end);

for p = 1:length(celldom)
    [idcelly idcellx] = find(masklabel == celldom(p));
    CoM(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
end

Nx = ACQinfo.pixelsPerLine;
Ny = ACQinfo.linesPerFrame;

trial = 3;
k = 1;
for p = 1:length(celldom)
    tc1 = cellS.cellMat{trial}(nID(p),:);
    tc1 = LFPfilt(tc1,0,1000/acqPeriod,8,.1);
    
    for q = p+1:length(celldom)
        
        tc2 = cellS.cellMat{trial}(nID(q),:);
        tc2 = LFPfilt(tc2,0,1000/acqPeriod,8,.2);
   
        dy = (CoM(p,1)-CoM(q,1))*150/Ny;
        dx = (CoM(p,2)-CoM(q,2))*150/Nx;
        
        if dy < 1000

            D(k) = sqrt(dy.^2 + dx.^2);

            dum = corrcoef(tc1,tc2);
            R(k) = dum(1,2);

            k = k+1;

        end
        
    end
end

figure,scatter(D,R,'.')

Nbins = 5;
ptsbin = floor(length(D(:))/Nbins);
[DS id] = sort(D(:));
RS = R(id);
for i = 1:Nbins
    idsamp = ((i-1)*ptsbin+1):i*ptsbin;
    if i == Nbins
        idsamp = ((i-1)*ptsbin+1):length(D(:));
    end
    sampD = DS(idsamp);
    sampR = RS(idsamp);
    
    muD(i) = trimmean(sampD,10);
    muR(i) = trimmean(sampR,10);
    sigD(i) = nanstd(sampD)/sqrt(length(sampD));
    sigR(i) = nanstd(sampR)/sqrt(length(sampD));
end
    
hold on
errorbar(muD,muR,sigR,'k')

xlabel('microns')