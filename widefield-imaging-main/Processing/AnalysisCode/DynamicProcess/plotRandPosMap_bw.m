function bwmap = plotRandPosMap_bw

%Ian Nauhaus

global kernelsIm G_RChandles maskS ACQinfo

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  

%Compute the orimap
dim = size(kernelsIm);
dim2 = size(kernelsIm{1,1,1}); 

t1 = 50; t2 = 450;
[dum id1] = min(abs(taudom-t1));
[dum id2] = min(abs(taudom-t2));

xdom = 0:(dim(2)-1);
xdom = xdom-mean(xdom);

dori = 180/dim(1);
oridom = 0:dori:(180-dori);

clear Tens
for bw = 1:dim(4)
    Tens{bw} = 0;
    for ori = 1:dim(1)
        for x = 1:dim(2)

            Tens{bw} = Tens{bw} + mean(kernelsIm{ori,x,1,bw}(:,:,id1:id2),3).^2;
        end
    end
    Tens{bw} = sqrt(Tens{bw});
end

h = fspecial('gaussian', size(Tens{1}(:,:,1)), 1);
h = abs(fft2(h));

Tens{1} = ifft2(fft2(Tens{1}).*h); %black is first index
Tens{2} = ifft2(fft2(Tens{2}).*h);

bwmap = (Tens{1}-Tens{2})./(Tens{1}+Tens{2});
%bwmap = log2(Tens{1}./Tens{2});
%id = find(isinf(bwmap));
%bwmap(id) = 0;
figure,imagesc(bwmap,[prctile(bwmap(:),.1) prctile(bwmap(:),99.9)])
axis image
colorbar
%figure,imagesc(bwmap,[])
title('off-on')


