NT = getnotrials;

dtau = 250;
taudom = 0:dtau:2000;  %This needs to have an element at 0

% slow movement correction
CHs = GetTrialData([1 1 0 0],1);
temp = mean(CHs{2}(:,:,2:end-2),3);

%use GUI to get the mask
global bwCell1 bwCell2
[bwCell1 bwCell2] = MakeCellMask(8,.8,3,6,mean(CHs{1},3));
%%
clear countall kernall
for T = 1:NT
    T
    CHs = GetTrialData([1 1 0 0],T);  %load second channel because it usually makes for better motion correction
    
    [Px_fast Py_fast] = getTrialMotion3(CHs{2});
    CHs{1} = makeGeoTrx(CHs{1},Px_fast,Py_fast);
    
    imdum = mean(CHs{2}(:,:,2:end-2),3);
    [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation for this trial
    for z = 1:length(CHs{1}(1,1,:))
        CHs{1}(:,:,z) = circshift(CHs{1}(:,:,z),[-mbest -nbest]); %transform
    end
    
    cond = getcondrep(T);
    [kern count] = getrevcorrkernel2(CHs{1},bwCell1,taudom,cond);    
    
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

for T = 1:NT
    CHs = GetTrialData([1 1 0 0],T);
    
%     [Px_fast Py_fast] = getTrialMotion3(CHs{2});
%     CHs{1} = makeGeoTrx(CHs{1},Px_fast,Py_fast);
    
    imdum = mean(CHs{1}(:,:,2:end-2),3);
    [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation for this trial
    for z = 1:length(CHs{1}(1,1,:))
        CHs{1}(:,:,z) = circshift(CHs{1}(:,:,z),[-mbest -nbest]); %transform
    end
    
    cond = getcondrep(T);
    [kern count] = getrevcorrkernel2(CHs{1},bwCell1,taudom,cond);    
    
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
