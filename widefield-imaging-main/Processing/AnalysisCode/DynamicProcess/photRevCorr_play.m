function [outshift] = photRevCorr_play(tau);

global ACQinfo

mPL = ACQinfo.msPerLine/1000/2;
acqPeriod = ACQinfo.linesPerFrame*mPL*1000  %ms per acquired frame
Npix = ACQinfo.linesPerFrame*ACQinfo.pixelsPerLine;

trialT = 60000;

%oriseq = stimseq(0:pepgetnoconditions-1);  %Get stimulus sequence

for i = 1:3
    oriseq{i} = 20*round((rand(1,120)*180 - 10)/20);
end
oridom = unique([oriseq{1} oriseq{2} oriseq{3}]);

% if oridom(end) == 999
%     oridom = oridom(1:end-1);
% end

%%%%%%%%%%%%%%%%%%%
sp = tau(2)-tau(1);
last = -acqPeriod-sign(rem(acqPeriod,sp))*sp;
dumtau = fliplr(tau(1):-sp:last);
tauBuild = [dumtau tau(2:end)];  %appended domain

out = cell(1,2);
out{1} = cell(1,length(oridom));
for i = 1:length(oridom)
    out{1}{i} = zeros(Npix,length(tauBuild));  %Initialize for accumulation
end
out{2} = out{1};
countmat = zeros(length(oridom),length(tauBuild));

for chan = 1:2  %Loop through both channels because of Matlab memory limits
    for k = 1:20

        pepsetcondition(k-1)

%         if chan == 1
%             CHs = GetTrialData(k,1,[1 0 0]);
%         elseif chan == 2
%             CHs = GetTrialData(k,1,[0 1 0]);
%         end

        CHs{1} = rand(Npix,ceil(trialT/acqPeriod));

        Tf = 1000/pepParam('refresh');  %Frame period in ms
        hper = pepgetparam('h_period');
        hper = hper(1);
        %hper = 1;
        Tupdate = Tf*hper;

        respDomain = (0:length(CHs{1}(1,:))-1)*acqPeriod;

        %oriseqdum = oriseq{k}(1:hper:end);
        N = round(trialT/(Tupdate));
        oriseqdum = 20*round((rand(1,N)*180 - 10)/20);
        %oriseqdum = oriseq{k};

        for ori = 1:length(oridom)

            id = find(oriseqdum == oridom(ori));

            stimes = (id-1)*Tupdate; %Stimulus times

            for i = 1:length(stimes)

                idx = find(respDomain > (stimes(i)+tauBuild(1)) & respDomain <= (stimes(i)+tauBuild(end)));
                if ~isempty(idx)
                    domidx = round((respDomain(idx)-(stimes(i)+tauBuild(1)))/sp) + 1; %domain indices within time domain
                    out{chan}{ori}(:,domidx) = out{chan}{ori}(:,domidx) + CHs{1}(:,idx);
                    countmat(ori,domidx) = countmat(ori,domidx) + 1;
                end

            end
        end
    end
    clear CHs

end
countmat = countmat/2  %It was doubled by the 2 channels
id = find(countmat(:) == 0);
mean(countmat(:))

for i = 1:length(oridom)
    for j = 1:length(tauBuild)
        out{1}{i}(:,j) = out{1}{i}(:,j)/countmat(i,j);  
        out{2}{i}(:,j) = out{2}{i}(:,j)/countmat(i,j); 
    end
end
    

delPix = 1000*mPL/ACQinfo.pixelsPerLine; %ms per pixel
iddel = (0:(Npix-1))*delPix;
iddel = round(iddel/sp);
iddel = [0 diff(iddel)];
iddel = find(iddel~=0);  %These are the indices to make a new time shift
L = length(tau);
LB = length(tauBuild);

%Shift each pixel in time

%First chunk
alp = LB-L+1;
ome = alp+L-1;
for i = 1:length(oridom)
    outshift{1}{i}(1:iddel(1)-1,:) = out{1}{i}(1:iddel(1)-1,alp:ome);
    outshift{2}{i}(1:iddel(1)-1,:) = out{2}{i}(1:iddel(1)-1,alp:ome);
end

%Middle chunks
for j = 1:length(iddel)-1
    %shift = round(delPix*iddel(j)/sp);
    shift = j;
    alp = (LB-L+1)-shift;
    ome = alp+L-1;
    for i = 1:length(oridom)
        outshift{1}{i}(iddel(j):iddel(j+1)-1,:) = out{1}{i}(iddel(j):iddel(j+1)-1,alp:ome);  
        outshift{2}{i}(iddel(j):iddel(j+1)-1,:) = out{2}{i}(iddel(j):iddel(j+1)-1,alp:ome); 
    end
end

%Last chunk
alp = (LB-L+1)-length(iddel);
ome = alp+L-1;
for i = 1:length(oridom)
    outshift{1}{i}(iddel(end):Npix,:) = out{1}{i}(iddel(end):Npix,alp:ome);
    outshift{2}{i}(iddel(end):Npix,:) = out{2}{i}(iddel(end):Npix,alp:ome);
end
