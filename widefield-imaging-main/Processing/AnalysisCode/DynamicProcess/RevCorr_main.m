NT = getnotrials;

dtau = 250;
taudom = 0:dtau:2000;  %This needs to have an element at 0

% slow movement correction
CH = GetTrialData([1 1 0 0],1);
[Px_fast Py_fast] = getTrialMotion3(CH{2});
CH{2} = makeGeoTrx(CH{2},Px_fast,Py_fast);
temp = mean(CH{2}(:,:,2:end-2),3);  %Template for trial-to-trial motion correction

%%
clear CHs CH
for T = 1:NT
    T
    CHs{T} = GetTrialData([1 1 0 0],T);  %load second channel because it usually makes for better motion correction
    
    [Px_fast Py_fast] = getTrialMotion3(CHs{T}{2});
    CHs{T}{1} = makeGeoTrx(CHs{T}{1},Px_fast,Py_fast);
    
    imdum = mean(CHs{T}{2}(:,:,2:end-2),3);
    [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation for this trial
    for z = 1:length(CHs{T}{1}(1,1,:))
        CHs{T}{1}(:,:,z) = circshift(CHs{T}{1}(:,:,z),[-mbest -nbest]); %transform
    end
    
end

%%


%use GUI to get the mask
global bwCell1 bwCell2

[bwCell1 bwCell2] = MakeCellMask(5,.5,3,5,mean(CHs{1}{1},3),3);

dtau = 250;
taudom = 0:dtau:2000;  %This needs to have an element at 0
kern = getrevcorrkernel(CHs,bwCell1,taudom);
kernelplots(kern,bwCell1,taudom);
%%

clear countall kernall 
for T = 1:length(CHs)
    T
     cond = getcondrep(T);
    [kern count] = getrevcorrkernel2(CHs{T}{1},bwCell1,taudom,cond);    
    
    if ~exist('countall')        
        for p = 1:length(kern)
            countall{p} = zeros(size(kern{p}));
            kernall{p} = zeros(size(kern{p}));
        end
    end
    
    for p = 1:length(kern)        
        countall{p} = countall{p} + count{p};      
        kernall{p} = kernall{p} + kern{p}; 
    end
    
end
%%
for p = 1:length(kern)
    
   kernall{p} = kernall{p}./countall{p}; 
    
end

kernelplots2(kernall,bwCell1,taudom);
