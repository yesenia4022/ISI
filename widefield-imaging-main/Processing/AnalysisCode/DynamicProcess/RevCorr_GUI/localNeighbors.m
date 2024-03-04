function localNeighbors

%Ian Nauhaus

% Local neighborhood analysis

global TC DM MK NB

NB = struct;  %make sure to reset

[xmicperpix ymicperpix] = getImResolution;


for c = 1:length(DM.colordom)
    %First normalize each tuning curve
    tcoriallmod = max(TC.tcoriall{c}')'*ones(1,length(TC.tcoriall{c}(1,:)));
    tcoriall = TC.tcoriall{c}./tcoriallmod;
    tcsfallmod = max(TC.tcsfall{c}')'*ones(1,length(TC.tcsfall{c}(1,:)));
    tcsfall = TC.tcsfall{c}./tcsfallmod;
    
    figure
    for p = 1:MK.Ncell

        dy = (MK.CoM(p,1) - MK.CoM(:,1))*ymicperpix;
        dx = (MK.CoM(p,2) - MK.CoM(:,2))*xmicperpix;
        r = sqrt(dx.^2+dy.^2); %all distances in microns        
               
        %%%Use distribution of preferences of local tuning curves to compute selectivity
        sig = 20; 
        wght = exp(-r.^2/(2*sig^2))'; 
        %wght = zeros(size(r')); id = find(r<sig); wght(id) = 1;
        id = find(~isnan(TC.sfpref{c}));
        sfmu = sum(wght(id).*log2(TC.sfpref{c}(id)))/sum(wght(id)); %first compute weighted mean of sf
        NB.sfvar_pks{c}(p) = sum(wght(id).*(log2(TC.sfpref{c}(id))-sfmu).^2)/sum(wght(id));  %weighted variance of log2(sfreq)
        NB.oriCV_pks{c}(p) = 1-abs(sum(wght(id).*exp(1i*2*TC.OAng{c}(id)*pi/180))/sum(wght(id)));  %Here, we use the distribution of preferred oris
        
        %%%Use the "population tuning curve" to compute selectivity        
        sig = 20; 
        wght = exp(-r.^2/(2*sig^2)); 
        %wght = zeros(size(r)); id = find(r<sig); wght(id) = 1; 
        wght(p) = 0;  %we exclude the center neuron to prevent bias
        wght = wght'/sum(wght);       
        Poptcori = nanmean(tcoriall.*(wght'*ones(1,length(DM.oridomI)))); 
        Poptcori = Poptcori/sum(Poptcori);
        [NB.oriCV_Pop{c}(p) NB.oripref_Pop{c}(p)] = orifind(Poptcori,DM.oridomI'); 
        NB.oriCV_Pop{c}(p) = 1-NB.oriCV_Pop{c}(p);
        
        Poptcsf = nanmean(tcsfall.*(wght'*ones(1,length(DM.sfdomI)))); 
        Poptcsf = phi(Poptcsf);
        Poptcsf = Poptcsf/sum(Poptcsf);        
        musf = sum(Poptcsf.*log2(DM.sfdomI));  %center of mass
        NB.sfvar_Pop{c}(p) = sum(Poptcsf.*log2(DM.sfdomI).^2) - musf^2;  %variance: E(x^2) - (E(x))^2    
        
%         [param ffit varacc ffitI domI] = DoGfit(Poptcsf,DM.sfdomI);
%         
%         [ma idma] = max(ffitI); [mi idmi] = min(ffitI(idma:end)); idmi = idmi+idma-1;
%         fpk = domI(idma);
%         thresh = (ma-mi)*.7 + mi;
%         [dum fhidum] = min(abs(ffitI(idma:idmi) - thresh));
%         fhi{c}(p) = domI(fhidum+idma-1);  %high cutoff sfreq
% 
%         [dum flodum] = min(abs(ffitI(1:idma) - thresh));
%         flo{c}(p) = domI(flodum);  %high cutoff sfreq
% 
%         NB.Qfac_Pop{c}(p) = fpk/(fhi{c}(p)-flo{c}(p));
%         NB.sfBW_Pop{c}(p) = log2(fhi{c}(p)/flo{c}(p));
%         NB.LPness_Pop{c}(p) = ffitI(1)/ffitI(idma);
%         
%         subplot(ceil(sqrt(MK.Ncell)),ceil(sqrt(MK.Ncell)),p)
%         semilogx(DM.sfdomI,Poptcsf), hold on
%         plot(domI,ffitI,'k')       
%         axis tight

    end
end


% figure,
% corrcoef(NB.oriCV_Pop{1},TC.OMag{1})
% subplot(2,2,1),scatter(NB.oriCV_Pop{1},TC.OMag{1},'.k'), xlabel('orimap circ var'), ylabel('ori selectivity')
% subplot(2,2,2), scatter(NB.sfvar_Pop{1},TC.sfBW{1},'.k'), xlabel('sfmap variance'), ylabel('sf selectivity')
% 
% subplot(2,2,3),scatter(NB.oriCV_Pop{1},NB.sfvar_Pop{1},'.k'), xlabel('orimap coherence'), ylabel('sfmap coherence')
% subplot(2,2,4), scatter(TC.OMag{1},TC.sfBW{1},'.k'), xlabel('ori selectivity'), ylabel('sf selectivity')

