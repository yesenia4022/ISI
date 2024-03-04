function [CHs] = getExptMean(Ntrials)

N = 0;

CHs = 0;
for trial = 1:Ntrials
    
    CHdum = GetTrialData([-inf inf],trial);
    
    N = N + length(CHdum(1,1,2:end-2));
    
    CHs = sum(CHdum(:,:,2:end-2),3) + CHs;
    
%     if chvec(1)
%         CHs{1} = sum(CHdum{1}(:,:,2:end-2),3) + CHs{1};
%     end
%     if chvec(2)
%         CHs{2} = sum(CHdum{2}(:,:,2:end-2),3) + CHs{2};
%     end
end


CHs = CHs/N;

    