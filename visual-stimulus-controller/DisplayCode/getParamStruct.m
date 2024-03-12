function PStruct = getParamStruct(OL_flag)

global Pstate Mstate Analyzer

if OL_flag %On-line flag; gets parameters from values in paramSelect GUI
    
    for i = 1:length(Pstate.param)
        eval(['PStruct.' Pstate.param{i}{1} '= Pstate.param{i}{3} ;'])
    end
    
else %gets values from loaded Analyzer file
    
    for i = 1:length(Analyzer.P)
        eval(['PStruct.' Analyzer.P.param{i}{1} '= Analyzer.P.param{i}{3} ;'])
    end
    
end


