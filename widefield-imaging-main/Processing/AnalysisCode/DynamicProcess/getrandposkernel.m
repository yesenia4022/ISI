function [kern kernblank countmat countmatblank] = getrandposkernel(cellMat,synctimes,bwCell1,Ntau)

global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

dtau = acqPeriod+50;
taudom = 0:dtau:dtau*Ntau
%taudom2 = linspace(0,1600,9)


hper = getparam('h_per');
%hper = 1; %in my stimulus code, the sequences have only one value for each
%presentation, so you can't set hper to one

expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt]
%load(['C:\2p_data\' Analyzer.M.anim '\log_files\' expt])
load(['F:\neurostuff\log_files\' expt])
rseeds = eval(Analyzer.L.param{1}{2})

xdom_dum = domains.xdom;

%Insert a place-holder for the blanks... the sequence will have index
%values that are one longer than the length of the spatial frequency
%domain, which are the blanks.

blankProb = 0

Tf = 1000/frate;  %Frame period in ms (frate obtained from log file) 
Tupdate = Tf*hper;

NT = getnotrials;

for t = 1:NT
    cond = getcondrep(t);
    s = Analyzer.loops.conds{cond}.val{1};

    eval(['bwS = rseed' num2str(s) '.bwseq;']);
    eval(['xS = rseed' num2str(s) '.xseq;']);
    eval(['oriS = rseed' num2str(s) '.oriseq;']);    
     
    bwseq{t} =  domains.bwdom(bwS);  
    xseq{t} =  xdom_dum(xS);     
    oriseq{t} =  domains.oridom(oriS);  
    
    %insert NaN for blanks
    if blankProb > 0
        idb = find(xS == length(domains.xdom)+1);
        
        bwseq{t}(idb) =  NaN;
        xseq{t}(idb) =  NaN;  %this is redundant
        oriseq{t}(idb) =  NaN;
    end
 
end


%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
xdom = domains.xdom;
bwdom = domains.bwdom;

Ncell = length(celldom);
NT = getnotrials;

figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(xdom),length(bwdom),length(taudom));
    kern{p} = zeros(length(oridom),length(xdom),length(bwdom),length(taudom));
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
    

        tcourse = zscore(tcourse);
        
        tdom = (0:length(tcourse)-1)*acqPeriod;
        %tdom_pix = tdom + tau_xy - ACQinfo.stimPredelay*1000;   %time domain of the pixel relative to onset of first stimulus
        tdom_pix = tdom + tau_xy;
        

        oriseqdum = oriseq{T}(1:end);
        xseqdum = xseq{T}(1:end);       
        bwseqdum = bwseq{T}(1:end);  


        for ori = 1:length(oridom)
            for x = 1:length(xdom)
                for bw = 1:length(bwdom)
                    
                    id = find(oriseqdum == oridom(ori) & xseqdum == xdom(x) & bwseqdum == bwdom(bw));
                    
                    %stimes = (id-1)*Tupdate; %Stimulus times
                    stimes = synctimes{T}(id)*1000;
                    
                    
                    for i = 1:length(stimes)
                        
                        idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
                        idx1 = idx1(1);
                        tpiece = idx1:idx1+length(taudom)-1;
                        
                        if tpiece(1)>0 & tpiece(end)<length(tcourse)
                            kern{p}(ori,x,bw,:) = squeeze(kern{p}(ori,x,bw,:)) + squeeze(tcourse(tpiece))';
                            countmat{p}(ori,x,bw,:) = countmat{p}(ori,x,bw,:) + 1;
                        end
                        
                    end
                end
            end
        end
        
        id = find(isnan(oriseqdum));
        stimes = synctimes{T}(id)*1000;
        
        for i = 1:length(stimes)
            
            idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
            idx1 = idx1(1);
            tpiece = idx1:idx1+length(taudom)-1;
            
            if tpiece(1)>0 & tpiece(end)<length(tcourse)
                kernblank{p} = kernblank{p} + tcourse(tpiece);
                countmatblank{p} = countmatblank{p} + 1;
            end
        end
            
        
    end


    kern{p} = kern{p}./countmat{p};
    %kernblank{p} = kernblank{p}./countmatblank{p};


    %normer = ones(length(kern{p}(:,1)),1)*mean(kern{p});
    %kern{p} = kern{p}-normer;
    
    %kernplot = squeeze(mean(kern{p}(:,1:end-2,:,:,:),3));
    %kernplot = squeeze(mean(kernplot,3));
    
    %kernblank = kern{p}(end,:);
    
    kernplot = squeeze(mean(kern{p},3));
    kernplot = squeeze(mean(kern{p}(3:8),3));
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    %kernplot = squeeze(nanmean(kernplot,2));
    %blankmat = ones(length(kernplot(:,1)),1)*kernblank{p}(:)';

    %kernplot = (kernplot-blankmat);
    %kernplot = [kernplot; kernblank{p}(:)'];
    
    imagesc(kernplot)
    
    drawnow
end
