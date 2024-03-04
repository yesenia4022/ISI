global G_RChandles G_handles

set(G_RChandles.LPflag,'value',1);
set(G_RChandles.HPflag,'value',1);
set(G_RChandles.LPWind,'value',1);
set(G_RChandles.HPWind,'value',1);
set(G_RChandles.Lwidth,'string',50);
set(G_RChandles.Hwidth,'string',5000);



%% Measure skewness

global maskS cellS

clear animS exptS
animS{1} = 'ab2'; 
animS{2} = 'ab2'; 
animS{3} = 'ab2';
animS{4} = 'ab3'; 
animS{5} = 'ab3'; 
animS{6} = 'aa9'; 
animS{7} = 'ab4';

animS{8} = 'ab2'; 
animS{9} = 'ab3'; 
animS{10} = 'ab3';

exptS{1} = 'u000_014'; %Small spot, but nice orthogonal linear zones
exptS{2} = 'u000_068'; %Linear ori and sf region; clearly orthogonal
exptS{3} = 'u001_004'; %fast linear zone, orthogonal sfreq map
exptS{4} = 'u002_014'; %isodomain, nice spat freq map
exptS{5} = 'u002_041'; %isodomain, awesome spat freq map (above/below 2_14)
exptS{6} = 'u001_027';  %pinwheel, low spat freq
exptS{7} = 'u001_003'; %ori fracture, low spat freq

exptS{8} = 'u000_087'; %small region. linear zone and orth sfreq
exptS{9} = 'u002_027'; %isodomain, sp freq linear zone
exptS{10} = 'u004_016';


Zdom = [-4:.01:10];
alltraceHist_control = 0;
alltraceHist = 0;
clear skewcell skewcell_control
for ex = 1:3
    ex
    anim = animS{ex};
    expt = exptS{ex};

    %load expt
    set(G_handles.loadana,'string',anim)
    set(G_handles.loadexp,'string',expt)
    Gsetdirectories

    hh = makeTemporalfilter;
    
    maskroot = 'C:\2ph_code\Beta\Masks\monkey\';
    maskpath = [maskroot anim '_' expt];
    load(maskpath,'maskS')
    
    nID = getNeuronMask;  %get the index values for the neurons
    masklabel = bwlabel(maskS.neuronmask,4);
    celldom = unique(masklabel);
    Ncell = length(nID);
    
    %%%Figure out which control neurons overlapped too much with other neurons
    xshift = 6; yshift = 7;    
    prcoverlap = getcontroloverlap(xshift,yshift);
    
    id = find(prcoverlap > 1);
    nID(id) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%skew of neurons%%%%%%%
    traceroot = 'e:\Beta\cellTraces\monkey\new2\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    load(tracepath,'cellS')
    
    Ntrial = length(cellS.cellMat);
    %Ntrial = 15;
    
    skewcell{ex} = 0;
    for i = 1:Ntrial  %loop each trial
        tracedum = cellS.cellMat{i}(nID,:); %pick out the neurons
        Nt = length(tracedum(1,:));
        %rng = round(Nt/10) : round(Nt*9/10);  
        rng = 1:Nt;        
        for j = 1:length(nID)  %smooth each cell
           tracedum(j,:) = ifft(hh.*fft(tracedum(j,:)));           
           tracedum(j,:) = (tracedum(j,:) - mean(tracedum(j,rng)))/std(tracedum(j,rng)); 
        end

        dum = tracedum(:,rng);
        alltraceHist = alltraceHist + hist(dum(:),Zdom)/Ntrial/length(dum(:));

        skewcell{ex} = skewcell{ex} + skewness(tracedum(:,rng)')/Ntrial;
    end    

    
    %%%%%skew of control%%%%%%%%
    traceroot = 'e:\Beta\cellTraces\monkey\new\control\';
    tracepath = [traceroot anim '_' expt '_cellS'];
    load(tracepath,'cellS')
    
    skewcell_control{ex} = 0;
    for i = 1:Ntrial %loop each trial
        tracedum = cellS.cellMat{i}(nID,:);  %pick out the neurons
        for j = 1:length(nID)  %smooth each cell
           tracedum(j,:) = ifft(hh.*fft(tracedum(j,:)));           
           tracedum(j,:) = (tracedum(j,:) - mean(tracedum(j,rng)))/std(tracedum(j,rng)); 
        end
        
        dum = tracedum(:,rng);
        alltraceHist_control = alltraceHist_control + hist(dum(:),Zdom)/Ntrial/length(dum(:));
        
        skewcell_control{ex} = skewcell_control{ex} + skewness(tracedum(:,rng)')/Ntrial;
    end

    %%%%%%%%%%%%%%
end
    
%%

skewAll = [];
skewCAll = [];
for ex = 1:2
    skewAll = [skewAll skewcell{ex}];
    skewCAll = [skewCAll skewcell_control{ex}];
end
id = find(skewCAll > 0 & skewAll > 0);

hdom = -.5:.15:1.7;
figure,subplot(3,2,1)
h = hist(skewAll(id),hdom); cumh = cumsum(h/length(id));
bar(hdom,h/length(id)), xlim([hdom(1) hdom(end)])
title(['med = ' num2str(median(skewAll(id))) ';  mu = ' num2str(mean(skewAll(id)))]), set(gca,'tickdir','out')
subplot(3,2,3)
h = hist(skewCAll(id),hdom); cumhC = cumsum(h/length(id));
bar(hdom,h/length(id)), xlim([hdom(1) hdom(end)])
title(['med = ' num2str(median(skewCAll(id))) ';  mu = ' num2str(mean(skewCAll(id)))]), set(gca,'tickdir','out')
subplot(3,2,5)
plot(hdom,cumh,'k')
hold on
plot(hdom,cumhC)
ylim([0 1.2]), xlim([-.5 1.4])
legend('Actual','Control')

subplot(3,2,2)
hdom = -10:.1:10;
Sratio = real(log10(skewAll(id)./skewCAll(id)));
h = hist(Sratio,hdom);
bar(hdom,h/sum(h)), xlim([-2 2]), set(gca,'tickdir','out')
title(['med = ' num2str(nanmedian(Sratio)) ';  geomean = ' num2str(10.^nanmean(real(Sratio)))])
set(gca,'tickdir','out')

[h p] = ttest2(skewAll(id),skewCAll(id))
[h p] = ttest(Sratio)

subplot(3,2,4), scatter(skewAll(id),skewCAll(id),'.k')
xlabel('skewness (actual)'), ylabel('skewness (control)')
hold on, plot([0 1],[0 1])
%%

alltraceHist = alltraceHist/sum(alltraceHist);
sqrt(sum(alltraceHist.*(Zdom.^2)))
figure, subplot(2,1,1)
bar(Zdom,alltraceHist)

alltraceHist_control = alltraceHist_control/sum(alltraceHist_control);
sqrt(sum(alltraceHist_control.*(Zdom.^2)))
subplot(2,1,2),bar(Zdom,alltraceHist_control)