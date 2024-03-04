function logmap = GprocessLog2(f0,bw,hh)

%f0dum is the cell array returned from fmeanimage.m


global Analyzer symbolInfo logTens G_handles

symID = symbolInfo.ID(1);
symID2 = symbolInfo.ID(2);

dom2index = get(G_handles.secCollapse,'value') - 1;  %Use zero to indicate the mean

val2 = [];
if dom2index>0
    dom2 = eval(Analyzer.L.param{symID2}{2});
    val2 = dom2(dom2index);
end


nc = length(Analyzer.loops.conds);

%if blank exists, it is always the last condition
bflag = 0;
if strcmp(Analyzer.loops.conds{nc}.symbol,'blank') 
    bflag = 1;
end

if ~isempty(hh)
    for k = 1:length(f0)
        %For log maps, we must filter before combining the images from each
        %condition because they are combined non-linearly.
        %if a filter exists, use it...
        
        id = find(isnan(f0{k}));
        f0{k}(id) = nanmedian(f0{k}(:));
        %f0{k} = rand(size(f0{k})); %Control
        f0{k} = ifft2(abs(fft2(hh)).*fft2(f0{k}));
    end
end

Tenscond = zeros(size(f0{1},1),size(f0{1},2),(nc-bflag));
for i = 1:(nc-bflag)
    logcond(i) = Analyzer.loops.conds{i}.val{symID}; %spatfreq for each condition
    dom2cond(i) = Analyzer.loops.conds{i}.val{symID2};
    Tenscond(:,:,i) = f0{i};
end

%blank subtraction
mi = min(Tenscond,[],3);
for i = 1:(nc-bflag)
    if bflag == 1 %if a blank exists
        Tenscond(:,:,i) = (Tenscond(:,:,i)-f0{end});
    else
        Tenscond(:,:,i) = Tenscond(:,:,i)-mi;
    end
end


logdom = unique(logcond);
logTens = zeros(size(f0{1},1),size(f0{1},2),length(logdom));
for i = 1:length(logdom)
    if isempty(val2)
        id = find(logcond == logdom(i));
    else
        id = find(logcond == logdom(i) & dom2cond == val2);  %probably used to establish an eye
    end
    logTens(:,:,i) = mean(Tenscond(:,:,id),3); %mean over all oris
    %logTens(:,:,i) = max(Tenscond(:,:,id),[],3);
end


%%Get maps from center of mass of tuning curve at each pixel

ma = max(logTens,[],3);  %Do this first so that it has dF/F units (for later)
mi = min(logTens,[],3);
%mapmag = (ma-mi)./(ma+mi);
mapmag = ma-mi;

logDomTens = zeros(size(logTens));
for i = 1:length(logdom)
    logDomTens(:,:,i) = log2(logdom(i));
end

[mi dum] = min(logTens,[],3);
for i = 1:length(logdom)
    logTens(:,:,i) = (logTens(:,:,i)-mi);
end

sumMat = sum(logTens,3);
for i = 1:length(logdom)
    logTens(:,:,i) = logTens(:,:,i)./sumMat;
end

mappref = sum(logTens.*logDomTens,3); %center of mass
mappref = 2.^mappref;  %put back in the right units

% [dum mappref] = max(logTens,[],3);
% mappref = logdom(mappref);


%%

logmap = mapmag + 1i*mappref;  %For convenience I make this one complex image
