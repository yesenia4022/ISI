function [tcPop CoM pardom] = GetCellTC(f0dum,cellmask)


global popState bsflag
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition

[kernPop popBlank CoM pardom] = GetCellKernels(f0dum,cellmask);

for i = 1:2
    pepsetcondition(i)
    if ~pepblank        
        Np = length(pepANA.listOfResults{i}.values);  %get number of looping parameters
        for z = 1:Np
            if strcmp(pepANA.listOfResults{i+1}.symbols(z),popState.funcSymbol)
                paramID = z;  %The parameter to analyze directly
            end
        end
    end
end


switch Np
    
    case 1
        
        tcPop = kernPop;

    case 2

        tcpopdum = squeeze(mean(kernPop,paramID));
        [dum idma] = max(tcpopdum);
        if paramID == 1
            for p = 1:Ncell
                tcPop(:,p) = squeeze(kernPop(:,idma(p),p));
            end
        else
            for p = 1:Ncell
                tcPop(:,p) = squeeze(kernPop(idma(p),:,p));
            end
        end

%     case 3
%         max(kernPop,[],);

end


%If functionality is orientation, then combine across opposite directions:
if strcmp(popState.funcSymbol,'ori')
    
    tcPopdum = (tcPop(1:end/2,:)+tcPop(end/2+1:end,:))/2;
    tcPop = tcPopdum;
    
    pardom{paramID} = pardom{paramID}(1:end/2);

end


if ~bsflag
    %popBlank = min(tcPop);
    blankmat = (popBlank'*ones(1,length(tcPop(:,1))))';
    tcPop = (tcPop-blankmat)./blankmat;  %deltaF/F
end

