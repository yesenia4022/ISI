function kern = photRevCorr2(CHs,shiftflag)

%2 does individual regions of interest

global ACQinfo

dtau = 50;
taudom = -3000:dtau:4800;

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

%oriseq = stimseq(0:pepgetnoconditions-1);  %Get stimulus sequence
oriseq = stimseq(0:4);  %Get stimulus sequence

oridom = unique([oriseq{1} oriseq{2} oriseq{3}]);

% if oridom(end) == 999
%     oridom = oridom(1:end-1);
% end

%%%%%%%%%%%%%%%%%%%

% for k = 1:length(oriseq)
%     k
%     %pepsetcondition(k-1)
% 
%     CHs{k} = GetTrialData([1 0 0],k);
%     
% end

%movement correction

% if shiftflag
%     temp = mean(CHs{1}{1}(:,:,10:30),3);
% 
%     for trial = 1:length(CHs)
%         for z = 1:length(CHs{trial}{1}(1,1,:))
% 
%             imdum = CHs{trial}{1}(:,:,z);  %Use first channel for alignment
% 
%             [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation
%             CHs{trial}{1}(:,:,z) = circshift(CHs{trial}{1}(:,:,z),[-mbest -nbest]); %transform
% 
%         end
%     end
% end

locs = [35 119; 73 90; 18 65; 83 47; 64 38; 103 40; 22 123; 51 56; 45 102];
Npos = length(locs(:,1));

W = 10;
figure
for p = 1:Npos
    xran = (locs(p,2)-floor(W/2)):(locs(p,2)+floor(W/2));
    yran = (locs(p,1)-floor(W/2)):(locs(p,1)+floor(W/2));

    tau_xy = (locs(p,1)-1)*ACQinfo.msPerLine + ptime*locs(p,2);

    countmat{p} = zeros(length(oridom),length(taudom));
    kern{p} = zeros(length(oridom),length(taudom));

    for k = 1:length(CHs)

        tcourse = squeeze(mean(mean(CHs{k}{1}(yran,xran,:),1),2));
        tcourse = LFPfilt(tcourse',0,1000/acqPeriod,inf,0.05)';

        %     hW = 21;
        %     hh = zeros(1,length(tcourse));
        %     hh(1:hW) = ones(1,hW);
        %     hh = -hh/sum(hh);
        %     hh(ceil(hW/2)) = hh(11) + 1;
        %     tcourse = ifft(fft(tcourse).*abs(fft(hh')));

        tcourse = zscore(tcourse);

        %Tf = 1000/pepParam('refresh');  %Frame period in ms (CRT)
        Tf = 1000/59.94;  %Frame period in ms  (LCD monitor)

        hper = pepgetparam('h_period');
        hper = hper(1);
        hper = 1;
        Tupdate = Tf*hper;

        tdom = (0:length(tcourse)-1)*acqPeriod;
        tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000 - 50;   %time domain of the pixel relative to onset of first stimulus

        oriseqdum = oriseq{k}(1:hper:end);

        for ori = 1:length(oridom)

            id = find(oriseqdum == oridom(ori));

            stimes = (id-1)*Tupdate; %Stimulus times

            for i = 1:length(stimes)
                
                idx = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<=stimes(i)+taudom(end)+dtau/2);
                
                domidx = round((tdom_pix(idx)-stimes(i))/dtau) + 1;
                domidx = domidx-taudom(1)/dtau;

                kern{p}(ori,domidx) = kern{p}(ori,domidx) + tcourse(idx)';
                countmat{p}(ori,domidx) = countmat{p}(ori,domidx) + 1;

            end
        end
    end
    kern{p} = kern{p}./countmat{p};
    
    subplot(ceil(sqrt(Npos)),ceil(sqrt(Npos)),p)
    imagesc(taudom,oridom(1:end-1),kern{p}(1:end-1,:))
    drawnow
end


for p = 1:Npos
    figure,
    subplot(2,2,1), imagesc(taudom,oridom(1:end-1),kern{p}(1:end-1,:))
    xlabel('time (ms)'), ylabel('ori')
    kernplot = kern{p}(1:end-1,:);

    [y x] = find(kernplot == max(kernplot(:)));

    tc1 = 0; tc2 = 0;
    tc = mean(kernplot(:,x+tc1:x+tc2),2);
    subplot(2,2,2), plot(taudom,kernplot(y,:)), hold on, plot(taudom,kernplot(y,:),'o'), xlabel('time (ms)')
    subplot(2,2,4), plot(oridom(1:end-1),tc), hold on, plot(oridom(1:end-1),tc,'o'), xlabel('ori')

end
