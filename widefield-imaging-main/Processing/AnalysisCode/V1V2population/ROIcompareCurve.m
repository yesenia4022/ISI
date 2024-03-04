function ROIcompareCurve

%compare Tuning curves from multiple ROIs

%ROIs and maps are cell arrays.  e.g. 1 cell array for V1 and 1 cell array for V2

global popState f0m bwCell1 pepANA

hh = [];

%%%%First get the tuning curves and blank response%%%%
[tcPopAll CoM pardom] = GetCellTC(f0m,bwCell1);
%%%%%%%%%%%%%%%%%%%%%

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

ROIs = popState.bwPopAnalysis;


NROI = length(ROIs);
for i = 1:NROI
     
    idRoi = inpolygon(CoM(:,2),CoM(:,1), popState.ROIPolyx{i}, popState.ROIPolyy{i});
    idRoi = find(idRoi);
  
    tcPop{i} = tcPopAll(:,idRoi)';  %Get subpopulation of tuning curves

    %Get rid of the "unresponsive" cells
    maxdF{i} = max(tcPop{i}');  %Get max dF/F for each cell
    id = find(maxdF{i} < popState.dFThresh);
    tcPop{i}(id,:) = [];
    NcellR = length(tcPop{i}(:,1));  %no. of responsive cells
         
       
    %Get population stats
    tcPopMod{i} = zeros(size(tcPop{i}));
    for p = 1:NcellR
        
        tcdum = tcPop{i}(p,:);
        
        if strcmp(popState.funcSymbol,'ori')
            Res = sum(tcdum.*exp(1i*2*pardom{paramID}*pi/180));
            Res = Res/sum(tcdum);
            ang = angle(Res);
            Opt{i}(p) = (ang+pi*(1-sign(ang)))/2*180/pi;
            Mag{i}(p) = abs(Res);
        else
            [ma Opt{i}(p)] = max(tcdum);
            Opt{i}(p) = pardom{paramID}(Opt{i}(p));
            Mag{i}(p) = ma;
        end

        if popState.alignflag
            idma = find(tcdum == max(tcdum));
            tcdum = circshift(tcdum,[0 round(length(pardom{paramID})/2)-idma]);            
        end
        
        if popState.peakflag
            tcdum = tcdum/max(tcdum);
        end
        
        tcPopMod{i}(p,:) = tcdum;
               
    end
end



colorid = 'brgky';

figure,
for i = 1:NROI
    muTC{i} = mean(tcPopMod{i});
    sigTC{i} = std(tcPopMod{i});
    
    errorbar(pardom{paramID},muTC{i},sigTC{i}/sqrt(length(tcPop{i}(:,1))),colorid(i))
    
    plotstr = popState.funcSymbol;
    plotstr(find(plotstr == '_')) = [];    
    xlabel(plotstr)
    if ~strcmp(popState.funcSymbol,'ori')
        set(gca,'Xscale','log')
    end
    hold on
end

% for i = 1:NROI
%     id = find(Mag{i}>1.5 | Mag{i}<0);
%     Mag{i}(id) = [];
%     Mag{i} = log(Mag{i});
% end

figure
for i = 1:NROI
    subplot(NROI,1,i)
    %dom = [0:.05:.4]+.025;
    [hp dom] = hist(maxdF{i},10);
    bar(dom,hp)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',colorid(i),'EdgeColor','w')
    xlabel('deltaF/F'), ylabel('N cells')
end

if strcmp(popState.funcSymbol,'ori')
    figure
    for i = 1:NROI
        subplot(NROI,1,i)
        dom = [0:.1:.9]+.05;
        hp = hist(Mag{i},dom);
        bar(dom,hp)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',colorid(i),'EdgeColor','w')
        xlabel('magnitude'), ylabel('N cells')
    end
end

figure
for i = 1:NROI
    subplot(NROI,1,i)
    %dom = [0:.1:.9]+.05;
    hp = hist(Opt{i},pardom{paramID});
    bar(pardom{paramID},hp)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',colorid(i),'EdgeColor','w')
    xlabel('Peak location'), ylabel('N cells')

    if ~strcmp(popState.funcSymbol,'ori')
        set(gca,'Xscale','log')
    end
end


if NROI == 2
    h = ttest2(Mag{2},Mag{1})

    if h == 1
        sprintf('Significantly different mean selectivity')
    else
        sprintf('Not significantly different mean preference')
    end
end


if NROI == 2
    h = ttest2(Opt{2},Opt{1})

    if h == 1
        sprintf('Significantly different mean preference')
    else
        sprintf('Not significantly different mean preference')
    end
end





