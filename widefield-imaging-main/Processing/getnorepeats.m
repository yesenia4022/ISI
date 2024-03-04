function nr = getnorepeats(cond)

global Analyzer

nr = length(Analyzer.loops.conds{cond}.repeats);