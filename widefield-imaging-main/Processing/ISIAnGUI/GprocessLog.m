function logmap = GprocessLog(f0,bw,hh,varargin)

%f0dum is the cell array returned from fmeanimage.m


global Analyzer symbolInfo logTens

symID = symbolInfo.ID(1);

if length(symbolInfo.ID) > 2
    symID2 = symbolInfo.ID(2);
end

nc = length(Analyzer.loops.conds);

%if blank exists, it is always the last condition
bflag = 0;
if strcmp(Analyzer.loops.conds{nc}.symbol,'blank') 
    bflag = 1;
end

%this is the id used to limit the tuning curves to some other value within
%the "second" parameter dimension (second dimension established in the GUI)
if ~isempty(varargin)
    secid = varargin{1};
else
    secid = [];
end

% for i = 1:length(Analyzer.loops.conds{1}.symbol)
%     if strcmp(Analyzer.loops.conds{1}.symbol{i},'s_freq')
%         symID = i;
%     end
% end

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
    Tenscond(:,:,i) = f0{i};
    
    if exist('symID2')
        seccond(i) = Analyzer.loops.conds{i}.val{symID2}; %Secondary parameter for each condition
    end
end

%blank subtraction
mi = min(Tenscond,[],3);
for i = 1:(nc-bflag)
    if bflag == 1 %if a blank exists
        %Tenscond(:,:,i) = (Tenscond(:,:,i)-f0{end});
    else
        Tenscond(:,:,i) = Tenscond(:,:,i)-mi;
    end
end

secmean = 1;

if secmean %mean over other parameters

    %Re-weight the sfreq tuning curve at each orientation by the average
    %response at that orientation
%     secdom = unique(seccond);
%     wt = zeros(size(f0{1},1),size(f0{1},2),length(secdom));
%     for i = 1:length(secdom)
%         id = find(seccond == secdom(i));
%         wt(:,:,i) = mean(Tenscond(:,:,id),3); %mean over all sfreqs
%     end
%     swt = sum(wt,3);
%     for i = 1:size(wt,3)
%        wt(:,:,i) = wt(:,:,i)./swt; %normalize the weighting function
%     end
%     for i = 1:length(secdom) %loop each ori
%         id = find(seccond == secdom(i));
%         for j = 1:length(id) %loop each sf index within this ori
%             Tenscond(:,:,id(j)) = Tenscond(:,:,id(j)).*wt(:,:,i); 
%         end
%     end


    logdom = unique(logcond);
    logTens = zeros(size(f0{1},1),size(f0{1},2),length(logdom));
    for i = 1:length(logdom)
        if isempty(secid)
            id = find(logcond == logdom(i));
        else
            id = find(logcond == logdom(i) & seccond == secid);  %probably used to establish an eye
        end
        logTens(:,:,i) = mean(Tenscond(:,:,id),3); %mean over all oris
        %logTens(:,:,i) = max(Tenscond(:,:,id),[],3);
    end

else %@ max
    
    secdom = unique(seccond);
    logpref = zeros(size(bw));
    logmag = zeros(size(bw));
    tclogall = 0;
    for i = 1:length(bw(:,1))
        
        for j = 1:length(bw(1,:))

            if bw(i,j)
                condlog = squeeze(Tenscond(i,j,:));
                [ma idsec] = max(condlog);                
                
                id = find(seccond == seccond(idsec));
                tclog = condlog(id);  %limit to the preferred secondary variable
                [logdom idsort] = sort(logcond(id));
                tclog0 = tclog(idsort);
                
                %also add 180
                idsec = idsec + length(secdom)/2;
                if idsec>length(secdom)
                    idsec = idsec-length(secdom);
                end
                id = find(seccond == seccond(idsec));
                tclog = condlog(id);  %limit to the preferred secondary variable
                [logdom idsort] = sort(logcond(id));
                tclog180 = tclog(idsort);                     
                
                tclog = (tclog0+tclog180)/2;
                
                tclog = tclog-min(tclog);
                tclog = tclog/sum(tclog);

                logTens(i,j,:) = tclog;

            end

        end
    end

end

%%
Gfitflag = 0;
if Gfitflag
    dlog = logdom(2)/logdom(1);
    mappref = zeros(size(bw));
    mapmag = zeros(size(bw));
    tclogall = 0;
    for i = 1:length(bw(:,1))
        i
        for j = 1:length(bw(1,:))

            if bw(i,j)
                tclog = squeeze(logTens(i,j,:));
                tclogall = tclogall+tclog;

                %             [param ffit varacc ffitI domI pk BW] = DoGfit(tclog',logdom);
                %             mappref(i,j) = pk;
                %             mapmag(i,j) = (max(ffitI)-min(ffitI))/(max(ffitI)+min(ffitI));

                [param ffit] = Gaussfit(log2([logdom dlog*logdom(end)]),[tclog' 0],0);

                mappref(i,j) = logdom(1)*dlog.^param(1); %this is how I get back to cyc/deg
                mapmag(i,j) = (max(ffit)-min(ffit))/(max(ffit)+min(ffit));
            end

        end
    end
    
else %Get maps from center of mass of tuning curve at each pixel
    
    ma = max(logTens,[],3);  %Do this first so that it has dF/F units (for later)
    mi = min(logTens,[],3);
    %mapmag = (ma-mi)./(ma+mi);
    mapmag = ma-mi;
    mapmag = ma;
    
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

    %[dum mappref] = max(logTens,[],3);
    %mappref = logdom(mappref);
    
end

%%

logmap = mapmag + 1i*mappref;  %For convenience I make this one complex image
