function [kern countmat] = getrevcorrkernel2(CHs,bwCell1,taudom,cond)


global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

dtau = taudom(2)-taudom(1);

Tf = 1000/100;  %Frame period in ms  (LCD monitor)

hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one
Tupdate = Tf*hper;

blankProb = getparam('blankProb');

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt]
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])

sfdom_dum = domains.sfdom;
%Insert a place-holder for the blanks... the sequence will have index
%values that are one longer than the length of the spatial frequency
%domain, which are the blanks.
if blankProb > 0
    sfdom_dum = [sfdom_dum NaN];   
end

rseed = Analyzer.loops.conds{cond}.val{1};

eval(['colorS = rseed' num2str(rseed) '.colorseq;']);
eval(['phaseS = rseed' num2str(rseed) '.phaseseq;']);
eval(['sfS = rseed' num2str(rseed) '.sfseq;']);
eval(['oriS = rseed' num2str(rseed) '.oriseq;']);

colorseq = domains.colordom(colorS);
phaseseq =  domains.phasedom(phaseS);
sfseq =  sfdom_dum(sfS);
oriseq =  domains.oridom(oriS);

%insert NaN for blanks
if blankProb > 0
    idb = find(sfS == length(domains.sfdom)+1);
    
    colorseq(idb) = NaN;
    phaseseq(idb) =  NaN;
    sfseq(idb) =  NaN;  %this is redundant
    oriseq(idb) =  NaN;
end


% if round(oridom(end)) == 999
%     oridom = oridom(1:end-1);
% end

%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;
colordom = domains.colordom;

Ncell = length(celldom);

for p = 1:Ncell
    
    [idcelly idcellx] = find(masklabel == celldom(p));
    idcell = find(masklabel(:) == celldom(p));
    
    CoM = [mean(idcelly) mean(idcellx)];  %center of mass
    
    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);
    
    countmat{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));
    kern{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));
    
    
    clear tcourse
    for z = 1:length(CHs(1,1,:))-1 %don't use last frame
        CHsdum = CHs(:,:,z);
        %tcourse(z) = (mean(CHsdum(idcell)) - mean(CHsdum(:)))/std(CHsdum(:));
        tcourse(z) = mean(CHsdum(idcell));
        
    end
    
    
    %         fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
    %         fdom = fdom(1:end-1);
    %         figure,plot(fdom,abs(fft(tcourse-mean(tcourse))))
    
    tcourse = LFPfilt(tcourse(:)',0,1000/acqPeriod,5,.05);
    %tcourse = BandErase(tcourse(:)',1000/acqPeriod,.66,.43);  %Breathing artifact
    
    %         hold on,plot(fdom,abs(fft(tcourse)),'r')
    %         asdf
    
    %     hW = 21;
    %     hh = zeros(1,length(tcourse));
    %     hh(1:hW) = ones(1,hW);
    %     hh = -hh/sum(hh);
    %     hh(ceil(hW/2)) = hh(11) + 1;
    %     tcourse = ifft(fft(tcourse).*abs(fft(hh')));
    
    tcourse = zscore(tcourse);
    
    tdom = (0:length(tcourse)-1)*acqPeriod;
    tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000 - 50;   %time domain of the pixel relative to onset of first stimulus

    
    for ori = 1:length(oridom)
        for sf = 1:length(sfdom)
            for phase = 1:length(phasedom)
                for color = 1:length(colordom)
                    
                    id = find(oriseq == oridom(ori) & sfseq == sfdom(sf) & phaseseq == phasedom(phase) & colorseq == colordom(color));
                    stimes = (id-1)*Tupdate; %Stimulus times
                    
                    for i = 1:length(stimes)
                        
                        for tauid = 1:length(taudom)
                            
                            idx = find(tdom_pix>=stimes(i)+taudom(tauid)-dtau/2 & tdom_pix<stimes(i)+taudom(tauid)+dtau/2);
                            
                            if ~isempty(idx)
                                kern{p}(ori,sf,phase,color,tauid) = kern{p}(ori,sf,phase,color,tauid) + mean(tcourse(idx));
                                countmat{p}(ori,sf,phase,color,tauid) = countmat{p}(ori,sf,phase,color,tauid) + 1;
                            end
                        end
                        
                    end
                end
            end
        end
    end
    

end
