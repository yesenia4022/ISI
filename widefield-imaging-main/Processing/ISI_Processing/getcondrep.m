function [cond rep] = getcondrep(ttag)

global Analyzer

conds = length(Analyzer.loops.conds);

for c = 1:conds
    reps = length(Analyzer.loops.conds{c}.repeats);
    for r = 1:reps
        ttag2 = Analyzer.loops.conds{c}.repeats{r}.trialno;
        if ttag == ttag2
            cond = c;
            rep = r;
            return
        end
    end
end

