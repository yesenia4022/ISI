function out = blankflag

global Analyzer

nc = getnoconditions;

if strcmp(Analyzer.loops.conds{nc}.symbol{1},'blank');
    out = 1;
else
    out = 0;
end