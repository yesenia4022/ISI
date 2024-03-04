function ROIcompareTime

%compare Time Courses from multiple ROIs

%ROIs and maps are cell arrays.  e.g. 1 cell array for V1 and 1 cell array for V2

global Tens popState bwCell1

hh = [];


NROI = length(ROIs);
Ncond = length(Tens);

ROIs = popState.bwPopAnalysis;

figure
for c = 1:Ncond
    for i = 1:NROI

        cmdum = ROIs{i}.*bwCell1;  %Cell mask within this ROI
        masklabel = bwlabel(cmdum);
        celldom = unique(masklabel);
        celldom = celldom(1:end);
        Ncell = length(celldom);

        %Get the time courses
        for p = 1:Ncell

            idcell = find(masklabel(:) == celldom(p));
            
%             [idcelly idcellx] = find(masklabel == celldom(p));
%             CoM{i,c}(p,:) = [mean(idcelly) mean(idcellx)];  %center of mass
            
            for tau = 1:length(Tens{c}(1,1,:))
                dum = Tens{c}(:,:,tau);
                tcPop{i,c}(p,tau) = mean(dum(idcell));
            end

        end
        
        subplot(Ncond,1,c)
        
        if i == 1
            plot(mean(tcPop{i,c}(:,1:end-1)))
        else
            plot(mean(tcPop{i,c}(:,1:end-1)),'r')
        end
        
        hold on
        
    end
end

