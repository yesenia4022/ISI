function [domains seqs] = getSeqInfo(trialdom)


global Analyzer G_RChandles

blankProb = getparam('blankProb');

logfileroot = get(G_RChandles.logfilePath,'string');
expt = [Analyzer.M.anim '_' Analyzer.M.unit '_' Analyzer.M.expt];

load([logfileroot Analyzer.M.anim '\' expt]')

if exist('rseed1')
    newSeq = 0;
elseif exist('randlog_T1')
    newSeq = 1;
end



%%%%For old sequence generation%%%%%%%
if ~newSeq  %get sequence and domains for old generation
    
    nparam = length(Analyzer.loops.conds{1}.symbol);
    rsbit = 0;
    for q = 1:length(nparam)
        if strcmp(Analyzer.loops.conds{1}.symbol{q},'rseed')
            rseedid = q;
            rseeds = eval(Analyzer.L.param{rseedid}{2});
            rsbit = 1;
        end
    end
    if ~rsbit
        rseeds = getparam('rseed');
        'Warning:  only 1 seed in this experiment'
    end
    
    sfdom_dum = domains.sfdom;
    %Insert a place-holder for the blanks... the sequence will have index
    %values that are one longer than the length of the spatial frequency
    %domain, which are the blanks.
    if blankProb > 0
        sfdom_dum = [sfdom_dum NaN];
    end
    
    for s = 1:length(rseeds)
        
        if exist(['rseed' num2str(s)])
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
        
    end    
    
    
    for q = 1:length(trialdom) %This ends up redundant, but it makes things consistent with the new sequence generation
        T = trialdom(q);
        [cond rep] = getcondrep(T);       
        
        if ~rsbit
            seedno = 1;
        else
            seedno = Analyzer.loops.conds{cond}.val{rseedid};
        end
        
        seqs{T}.oriseq = oriseq{seedno};  
        seqs{T}.sfseq = sfseq{seedno};
        seqs{T}.phaseseq = phaseseq{seedno};
        seqs{T}.colorseq = colorseq{seedno};
       
        domainsdum{T} = domains;
        
    end
    domains = domainsdum;
    clear domainsdum
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%For new sequence generation%%%%%%%

if newSeq  
    
    for q = 1:length(trialdom)      
        
        T = trialdom(q);
        
        eval(['seq = randlog_T' num2str(T) '.seqs;']);
        eval(['domains{' num2str(T) '} = randlog_T' num2str(T) '.domains;']);  %in the new versions the domain can change for each trial (e.g. via a looper variable)
        
        oriseq = domains{T}.oridom(seq.oriseq);
        phaseseq = domains{T}.phasedom(seq.phaseseq);
        colorseq = domains{T}.colordom(seq.colorseq);
        
        if blankProb > 0
            idb = find(seq.sfseq == length(domains{T}.sfdom)+1); %Blanks are identified with spatial frequency sequence
            oriseq(idb) = NaN;
            phaseseq(idb) = NaN;
            colorseq(idb) = NaN;
            
            sfdom_dum = [domains{T}.sfdom NaN];
            sfseqdum = sfdom_dum(seq.sfseq);
        else
            sfseqdum = domains{T}.sfdom(seq.sfseq);
        end
        
        seqs{T}.oriseq = oriseq;
        seqs{T}.sfseq = sfseqdum;
        seqs{T}.phaseseq = phaseseq;
        seqs{T}.colorseq = colorseq;
        
        
    end
    
end