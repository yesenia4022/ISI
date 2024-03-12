function pval = getParamVal(psymbol,OL_flag)

global Analyzer Pstate

if OL_flag %On-line flag; gets parameters from values in paramSelect GUI
    
    for i = 1:length(Pstate.param)
        if strcmp(psymbol,Pstate.param{i}{1})
            idx = i;
            break;
        end
    end
    
    pval = Pstate.param{idx}{3}
    
else %off-line; get it from loaded analyzer file
    
    for i = 1:length(Analyzer.P.param)
        if strcmp(Analyzer.P.param{i}{1},psymbol)
            
            pval = Analyzer.P.param{i}{3};
            break
            
        end
    end
    
end
