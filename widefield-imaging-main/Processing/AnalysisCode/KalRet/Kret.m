function Kret

global ACQinfo

stimfreq = pepgetparam('t_period');
stimfreq = 60/stimfreq(1);   %Hz

for c = 1:pepgetnoconditions
    pepsetcondition(c-1)
    for r = 1:pepgetnorepeats
        pepsetrepeat(r-1)
        CHs = GetTrialData([1 0 0]);
        
        for i = 1:length(CHs{1}(1,1,:))
            CHsdum = CHs{2}{1}(:,:,i)'; 
            CHsvec(:,i) = CHsdum(:)
        end
        clear CHs
        
        

        
    end
end

pos = [137 113];

W = 30;
xran = (pos(2)-floor(W/2)):(pos(2)+floor(W/2));
yran = (pos(1)-floor(W/2)):(pos(1)+floor(W/2));

global ACQinfo

for i = 1:length(CHs{2}{1}(1,1,:))
    CHsdum = CHs{2}{1}(yran,xran,i);
    tcourse(i) = mean(CHsdum(:));
end

sp = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
tdom = (0:length(tcourse)-1)*sp;
fdom = linspace(0,1,length(tcourse)+1)*1000/sp;
fdom = fdom(1:end-1);

tcourse = zscore(tcourse);
figure,plot(tdom/1000,tcourse)
figure,plot(fdom,abs(fft(tcourse)))


stimfreq = pepgetparam('t_period');
stimfreq = 60/stimfreq(1);   %Hz

CHsF = fft(CHs{2}{1},256,3);
fdom = linspace(0,1,length(CHsF(1,1,:))+1)*1000/sp;
[dum id] = min(abs(stimfreq-fdom));

ImageF1 = CHsF(:,:,id);

figure,imagesc(abs(ImageF1))
figure,imagesc(angle(ImageF1))
