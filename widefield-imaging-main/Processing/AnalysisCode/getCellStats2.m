function getCellStats2

global ACQinfo cellS bsflag G_handles

bsflag = get(G_handles.basesub,'value');  %Frame start in ms (to average)

Flim = str2double(get(G_handles.epistart,'String'));  %Frame start in ms (to average)
Flim(2) = str2double(get(G_handles.epistop,'String')); %Frame stop in ms (to average)

framePer = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %frame period in ms
Flim = Flim+ACQinfo.stimPredelay*1000;  %user input is relative to stimulus onset, not trial beginning

Frame1 = floor(Flim(1)/framePer) + 1;
Frame2 = ceil(Flim(2)/framePer) + 1;

Nt = length(cellS.cellMat{1}(1,:,1));
Ncell = length(cellS.cellMat{1}(:,1,1));
dom = linspace(-framePer*Nt/2,framePer*Nt/2,Nt);
sig = 300;  %ms
smoother = exp(-dom.^2/(2*sig^2));
smoother = smoother/sum(smoother);
smootherblk = ones(Ncell,1)*smoother;
smootherblk = abs(fft(smootherblk,[],2));

nf = length(Frame1:Frame2);

%Baseline normalization
if bsflag
    blim = str2double(get(G_handles.bstart,'String')); %in msec as well
    blim(2) = str2double(get(G_handles.bstop,'String'));
    blim = blim+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning
    bframe1 = floor(blim(1)/framePer) + 1;
    bframe2 = ceil(blim(2)/framePer) + 1;
    for i = 1:length(cellS.cellMat)  %loop through each condition
        cellS.cellMat_norm{i} = zeros(size(cellS.cellMat{i})); %Preallocate
        for j = 1:length(cellS.cellMat{i}(:,1,1)) %loop through each cell
            condcellMat = squeeze(cellS.cellMat{i}(j,:,:)); %[time repeats]
            blank = mean(condcellMat(bframe1:bframe2,:));
            blank = ones(length(condcellMat(:,1)),1)*blank;
            condcellMat = (condcellMat-blank)./blank;
            cellS.cellMat_norm{i}(j,:,:) = condcellMat;
        end
    end
    
end

%Stats
for i = 1:length(cellS.cellMat)  %loop through each condition
    
    nr = getnorepeats(i);
    
    dim = size(cellS.cellMat{i});
    
    if bsflag
        cellS.muTime{i} = squeeze(mean(cellS.cellMat_norm{i},3)); %mean across repeats
    else
        cellS.muTime{i} = squeeze(mean(cellS.cellMat{i},3)); %mean across repeats
    end
    cellS.mu{i} = mean(cellS.muTime{i}(:,Frame1:Frame2),2);  %mean across repeats and time window
    
    muTimesmooth = ifft(fft(cellS.muTime{i},[],2).*smootherblk,[],2);  %smooth before taking max
    cellS.maxi{i} = max(muTimesmooth(:,Frame1:Frame2),[],2);  %mean across repeats; max over time window

%     dum = reshape(cellS.cellMat{i}(:,Frame1:Frame2,:),dim(1),dim(3)*length(Frame1:Frame2));
%     cellS.sig{i} = std(dum,[],2)/sqrt(nr*nf); %std error across time and repeats
    
    %I think computing the mean across time first (as done below) is a
    %better way to compute the standard error for each condition.
    %Otherwise, the standard error becomes really low due to the number of
    %time samples (as above). It also doesn't really make sense to compute
    %a standard deviation across two dimensions (i.e. repeats and time
    %points)
    if bsflag
        cellS.sigTime{i} = std(cellS.cellMat_norm{i},[],3)/sqrt(nr); %std error across repeats
        dum = squeeze(mean(cellS.cellMat_norm{i}(:,Frame1:Frame2,:),2));
    else
        cellS.sigTime{i} = std(cellS.cellMat{i},[],3)/sqrt(nr); %std error across repeats
        dum = squeeze(mean(cellS.cellMat{i}(:,Frame1:Frame2,:),2));
    end
    cellS.sig{i} = std(dum,[],2)/sqrt(nr); %std error across repeats, for each condition
    
end



mu2 = cellS.mu;
%Compute stats
for p = 1:length(cellS.muTime{i}(:,1))
    
    clear tc tctime
    for i = 1:length(cellS.cellMat) 
        tc(i) = mu2{i}(p);
        tctime(i,:) = cellS.muTime{i}(p,:);
    end
    
    tctime = mean(tctime);
    tctime = ifft(fft(tctime.*smoother));

    [ma idpk] = max(tctime);
    
    %FRONT END
    tcdum = tctime(1:idpk);

    %time to 50%
    thresh = (ma+tctime(1))*.5;
    [dum id] = min(abs(tcdum-thresh));
    tt50 = id;


    %BACK END
    tcdum = tctime(idpk:end);

    %time to 50%
    thresh = (ma+tctime(1))*.7;
    [dum id] = min(abs(tcdum-thresh));
    id = id + idpk - 1;
    tt50back = id;
    
    for i = 1:length(cellS.cellMat)
        cellS.mu{i}(p) = mean(cellS.muTime{i}(p,tt50:tt50back),2);
    end
        

end