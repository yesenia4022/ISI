function out = stimblank(cond)

global Analyzer


if strcmp(Analyzer.loops.conds{cond}.symbol{1},'blank');
    out = 1;
else
    out = 0;
end