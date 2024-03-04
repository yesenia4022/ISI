function AnalyzerPatch

global Analyzer

nc = length(Analyzer.loops.conds);


if strcmp(Analyzer.loops.conds{nc}.symbol{1},'blank')
    
    for c = 1:nc-1  %we don't want ttagMat to contain the blanks
        nr = length(Analyzer.loops.conds{c}.repeats);
        for r = 1:nr
            ttagMat(c,r) = Analyzer.loops.conds{c}.repeats{r}.trialno;
        end
    end  

    Nbreps = length(Analyzer.loops.conds{nc}.repeats);
    %I know its wrong if the blank trials also exist in the normal
    %trials:
    if ~isempty(find(Analyzer.loops.conds{nc}.repeats{Nbreps}.trialno == ttagMat))
        sprintf('Assigned trial numbers for blanks is wrong in analyzer file. Applying patch (I.N.)...')
        
        blankcounter = 0;
        for i = 1:Nbreps  
            
            Analyzer.loops.conds{nc}.repeats{i}.trialno = Analyzer.loops.conds{nc}.repeats{i}.trialno + blankcounter;
            blankcounter = blankcounter+1;
            
        end
        
    end
    
end
