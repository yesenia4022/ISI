function countmat = getGratingPresentationCounts(trialdom)

[domains seqs] = getSeqInfo(trialdom);

oridom = domains{trialdom(1)}.oridom;
sfdom = domains{trialdom(1)}.sfdom;
phasedom = domains{trialdom(1)}.phasedom;
colordom = domains{trialdom(1)}.colordom;

countmat = zeros(length(oridom),length(sfdom),length(phasedom),length(colordom));
for trialid = 1:length(trialdom)

    T = trialdom(trialid);

    for ori = 1:length(oridom)
        for sf = 1:length(sfdom)
            for phase = 1:length(phasedom)
                for color = 1:length(colordom)

                    id = find(seqs{T}.oriseq == oridom(ori) & seqs{T}.sfseq == sfdom(sf) & seqs{T}.phaseseq == phasedom(phase)& seqs{T}.colorseq == colordom(color));
                    
                    countmat(ori,sf,phase,color) = countmat(ori,sf,phase,color) + length(id);

                end
            end
        end
    end
end


