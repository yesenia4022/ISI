function [tcmat tcourseHi tcourseLo primDom blank legStr] = getcelldata(pos,W)

%A modified get pixel data... this uses the cell mask.

global Analyzer symbolInfo cellS maskS

masklabel = bwlabel(maskS.bwCell1);
celldom = unique(masklabel);

cellID = masklabel(pos(2),pos(1));

nc = getnoconditions;

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
Nloop = nc;
blank = [];
if bflag
    Nloop = nc-1;    
    dum = f0m{end}(yran,xran);
    blank = mean(dum(:));
end

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping parameters

for i = 1:Nsym
    allDom{i} = getdomain(symbolInfo.str{i});
    dim(i) = length(allDom{i});
end

figure
for i = 1:length(cellS.muTime)   
    for j = 1:Nsym
        val = Analyzer.loops.conds{i}.val{j};
        idval(j) = find(allDom{j} == val);
    end    
    subplot(dim(1),dim(2),i)
    plot(cellS.muTime{i}(cellID,:))
end
