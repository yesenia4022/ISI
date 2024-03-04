function hper = gethper(cond)

global Analyzer

Nsym = length(Analyzer.loops.conds{cond}.symbol);

hper = getparam('h_per');  %Get the default
for i = 1:Nsym  %...or the one for this condition

    if strcmp('h_per',Analyzer.loops.conds{cond}.symbol{i});
        hper = Analyzer.loops.conds{cond}.val{i};
        break
    end
    
end

