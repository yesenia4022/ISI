function trial = gettrial(cond,rep)

global Analyzer

trial = Analyzer.loops.conds{cond}.repeats{rep}.trialno;
