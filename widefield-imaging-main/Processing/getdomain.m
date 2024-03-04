function domain = getdomain(symb)

global Analyzer

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping params
for idsym = 1:Nsym
    if strcmp(symb,Analyzer.loops.conds{1}.symbol{idsym});
        break
    end
end

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
Nloop = getnoconditions;
if bflag
    Nloop = Nloop-1;
end

for i = 1:Nloop
    domain(i) = Analyzer.loops.conds{i}.val{idsym};
end

domain = unique(domain);
