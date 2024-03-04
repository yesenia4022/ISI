clear CHs
pepsetcondition(0)
for r = 1:pepgetnorepeats
    pepsetrepeat(r-1)
    CHs{r} = GetTrialData([1 0 1]);
end



%movement correction
shiftflag = 1;
if shiftflag
    temp = mean(CHs{1}{1}(:,:,10:30),3);

    for trial = 1:length(CHs)
        for z = 1:length(CHs{trial}{1}(1,1,:))

            imdum = CHs{trial}{1}(:,:,z);  %Use first channel for alignment

            [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation
            CHs{trial}{1}(:,:,z) = circshift(CHs{trial}{1}(:,:,z),[-mbest -nbest]); %transform

        end
    end
end



bw = ones(size(CHs{1}{1}(:,:,1)));
[bwCell1 bwCell2] = MakeCellMask(14,.8,3);
bwCell1 = bwCell1.*bw;

[kern popResp CoM] = flashRep2(CHs,bwCell1);

[cc D] = flashPairwiseAnal(popResp,CoM);

plotccProfile(cc,bwCell1)




load('C:\Documents and Settings\Acquire\Desktop\flashRepeat\kern')
playTensor(kern,.1)



load('C:\Documents and Settings\Acquire\Desktop\flashRepeat\muStack')
playRedGreenTens(muStack,.2)

