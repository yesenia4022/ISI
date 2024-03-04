function [kern kernblank] = getrevcorrkernel4(cellMat,synctimes,bwCell1,taudom)


%3 takes the cell time courses as input, instead of all the images

%4 collapses over phase and color

global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

dtau = taudom(2)-taudom(1);



hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

blankProb = getparam('blankProb');

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt]
load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
%load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2});

sfdom_dum = domains.sfdom;
%Insert a place-holder for the blanks... the sequence will have index
%values that are one longer than the length of the spatial frequency
%domain, which are the blanks.
if blankProb > 0
    sfdom_dum = [sfdom_dum NaN];   
end


Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 
Tupdate = Tf*hper;

for i = 1:33
    
    eval(['colorS = rseed' num2str(i) '.colorseq;']);
    eval(['phaseS = rseed' num2str(i) '.phaseseq;']);
    eval(['sfS = rseed' num2str(i) '.sfseq;']);
    eval(['oriS = rseed' num2str(i) '.oriseq;']);    
    
    colorseq{i} = domains.colordom(colorS);   
    phaseseq{i} =  domains.phasedom(phaseS);  
    sfseq{i} =  sfdom_dum(sfS);     
    oriseq{i} =  domains.oridom(oriS);  
    
    %insert NaN for blanks
    if blankProb > 0
        idb = find(sfS == length(domains.sfdom)+1);
        
        colorseq{i}(idb) = NaN;
        phaseseq{i}(idb) =  NaN;
        sfseq{i}(idb) =  NaN;  %this is redundant
        oriseq{i}(idb) =  NaN;
    end
 
end


%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;
colordom = domains.colordom;

Ncell = length(celldom);
NT = getnotrials;

figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(sfdom),length(taudom));
    kern{p} = zeros(length(oridom),length(sfdom),length(taudom));
    countmatblank{p} = zeros(1,length(taudom));
    kernblank{p} = zeros(1,length(taudom));

    for T = 1:length(cellMat)

        tcourse = cellMat{T}(p,:);

        %         fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
        %         fdom = fdom(1:end-1);
        %         figure,plot(fdom,abs(fft(tcourse-mean(tcourse))))

        %mu = mean(tcourse);
        tcourse = LFPfilt(tcourse(:)',0,1000/acqPeriod,2,.05);
        %tcourse = tcourse+mu;  %replace the mean
     

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
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;
        
        %idbase = find(tdom_pix<synctimes{T}(1)*1000 | tdom_pix>(synctimes{T}(end)*1000+500));
        %bLine = mean(tcourse(idbase));
        
        %tcourse = (tcourse-bLine)/bLine;
        
        [cond] = getcondrep(T);

        oriseqdum = oriseq{cond}(1:end);
        sfseqdum = sfseq{cond}(1:end);       

        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)
                
                id = find(oriseqdum == oridom(ori) & sfseqdum == sfdom(sf));
                
                %stimes = (id-1)*Tupdate; %Stimulus times
                stimes = synctimes{T}(id)*1000;
                
                
                for i = 1:length(stimes)
                    
                    for tauid = 1:length(taudom)
                        stimes+taudom(tauid)-dtau/2
                        
                        idx = find(tdom_pix>stimes+taudom(tauid)-dtau/2 & tdom_pix<stimes+taudom(tauid)+dtau/2);
                        
                        if ~isempty(idx)
                            kern{p}(ori,sf,tauid) = kern{p}(ori,sf,tauid) + mean(tcourse(idx));
                            countmat{p}(ori,sf,tauid) = countmat{p}(ori,sf,tauid) + 1;
                        end
                    end
                end
            end
        end
        
        
        id = find(isnan(oriseqdum));
        stimes = synctimes{T}(id)*1000;
        
        for i = 1:length(stimes)
            
            for tauid = 1:length(taudom)
                
                idx = find(tdom_pix>stimes(i)+taudom(tauid)-dtau/2 & tdom_pix<stimes(i)+taudom(tauid)+dtau/2);
                
                if ~isempty(idx)
                    kernblank{p}(tauid) = kernblank{p}(tauid) + mean(tcourse(idx));
                    countmatblank{p}(tauid) = countmatblank{p}(tauid) + 1;
                end
            end
            
        end
        
    end

    
    kern{p} = kern{p}./countmat{p};
    kernblank{p} = kernblank{p}./countmatblank{p};

    
    kernplot = squeeze(mean(kern{p}(:,1:end-2,:),2));
    

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    imagesc(kernplot)
    %plot(kernplot(1:4,:)')
    
    
    drawnow
end
