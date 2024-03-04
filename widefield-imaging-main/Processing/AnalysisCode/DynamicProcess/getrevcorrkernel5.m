function [kern kernblank countmat countmatblank] = getrevcorrkernel5(cellMat,synctimes,bwCell1,Ntau)

%3 takes the cell time courses as input, instead of all the images

global ACQinfo Analyzer

%%%%

masklabel = bwlabel(bwCell1);
celldom = unique(masklabel);
celldom = celldom(1:end);

%%%%

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %ms per acquired frame
ptime = ACQinfo.msPerLine/ACQinfo.pixelsPerLine;  %pixel time (ms)

dtau = acqPeriod+50;
taudom = 0:dtau:dtau*Ntau;
%taudom2 = linspace(0,1600,9)


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


for s = 1:length(rseeds)
    
    eval(['colorS = rseed' num2str(s) '.colorseq;']);
    eval(['phaseS = rseed' num2str(s) '.phaseseq;']);
    eval(['sfS = rseed' num2str(s) '.sfseq;']);
    eval(['oriS = rseed' num2str(s) '.oriseq;']);    
    
    colorseq{s} = domains.colordom(colorS);   
    phaseseq{s} =  domains.phasedom(phaseS);  
    sfseq{s} =  sfdom_dum(sfS);     
    oriseq{s} =  domains.oridom(oriS);  
    
    %insert NaN for blanks
    if blankProb > 0
        idb = find(sfS == length(domains.sfdom)+1);
        
        colorseq{s}(idb) = NaN;
        phaseseq{s}(idb) =  NaN;
        sfseq{s}(idb) =  NaN;  %this is redundant
        oriseq{s}(idb) =  NaN;
    end
 
end


%%%%%%%%%%%%%%%%%%%
oridom = domains.oridom;
sfdom = domains.sfdom;
phasedom = domains.phasedom;
colordom = domains.colordom;

Ncell = length(celldom);
NT = getnotrials;
tcoursedum = 0;
figure
for p = 1:Ncell

    [idcelly idcellx] = find(masklabel == celldom(p));

    CoM = [mean(idcelly) mean(idcellx)];  %center of mass

    tau_xy = (CoM(1)-1)*ACQinfo.msPerLine + ptime*CoM(2);

    countmat{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));
    kern{p} = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom),length(taudom));
    countmatblank{p} = zeros(1,length(taudom));
    kernblank{p} = zeros(1,length(taudom));

    for T = 1:2:length(cellMat)

        tcourse = cellMat{T}(p,:);

         tcoursedum = abs(fft(mean(cellMat{T}))) + tcoursedum;
%                 fdom = linspace(0,1000/acqPeriod,length(tcourse)+1);
%                 fdom = fdom(1:end-1);
%                 figure,plot(fdom,abs(fft(tcoursedum-mean(tcoursedum))))

        %mu = mean(tcourse);
        tcourse = LFPfilt(tcourse(:)',0,1000/acqPeriod,4,.05);
        %tcourse = tcourse+mu;  %replace the mea
     

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
        
        seedno = Analyzer.loops.conds{cond}.val{1};

        oriseqdum = oriseq{seedno};
        sfseqdum = sfseq{seedno};       
        phaseseqdum = phaseseq{seedno};  
        colorseqdum = colorseq{seedno}; 
        
        

        for ori = 1:length(oridom)
            for sf = 1:length(sfdom)
                for phase = 1:length(phasedom)
                    for color = 1:length(colordom)
                        
                        id = find(oriseqdum == oridom(ori) & sfseqdum == sfdom(sf) & phaseseqdum == phasedom(phase)& colorseqdum == colordom(color));
                        
                        %stimes = (id-1)*Tupdate; %Stimulus times
                        stimes = synctimes{T}(id)*1000;
                        
                        
                        for i = 1:length(stimes)

                            idx1 = find(tdom_pix>=stimes(i)+taudom(1)-dtau/2 & tdom_pix<stimes(i)+taudom(1)+dtau/2);
                            idx1 = idx1(1);
                            tpiece = idx1:idx1+length(taudom)-1;
                            
                            if tpiece(1)>0 & tpiece(end)<length(tcourse)
                                kern{p}(ori,sf,phase,color,:) = squeeze(kern{p}(ori,sf,phase,color,:)) + squeeze(tcourse(tpiece))';
                                countmat{p}(ori,sf,phase,color,:) = countmat{p}(ori,sf,phase,color,:) + 1;
                            end
                            
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
    kernblank{p} = kernblank{p}./countmatblank{p};


    %normer = ones(length(kern{p}(:,1)),1)*mean(kern{p});
    %kern{p} = kern{p}-normer;

    kernplot = mean(kern{p}(:,1:end-1,:,:,:),3);  %mean across phase
    kernplot = mean(kernplot,4); %mean across color
    
    %kernblank = kern{p}(end,:);
    
    %kernplot = ifft(fft(kernplot).*abs(fft(kernsmooth)));

    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    
    kernplot = squeeze(nanmean(kernplot,2));  %mean across spatial freqency
    blankmat = ones(length(kernplot(:,1)),1)*kernblank{p}(:)';

    %kernplot = (kernplot-blankmat);
    %kernplot = [kernplot; kernblank{p}(:)'];
    
    imagesc(kernplot)
    
    drawnow
end


