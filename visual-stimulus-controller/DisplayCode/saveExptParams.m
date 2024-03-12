function saveExptParams

global Mstate Pstate Lstate looperInfo ACQ

%Save the analzer file

Analyzer.M = Mstate;
Analyzer.P = Pstate;
Analyzer.L = Lstate;
Analyzer.loops = looperInfo;

if Mstate.WF
    Analyzer.ACQ = ACQ;  %This is a new addition.  Should be saving this.
end

title = [Mstate.anim '_' sprintf('u%s',Mstate.unit) '_' Mstate.expt];

roots = parseString(Mstate.analyzerRoot,';');

%Also save on the network to 2p aquisition computer, if needed
if Mstate.twoP
    roots{length(roots)+1} = Mstate.analyzerRoot_2p;
end

%Save each root:
for i = 1:length(roots)

    dd = [roots{i} '\' Mstate.anim];
dd
    if(~exist(dd))
        mkdir(dd);  %if there is a new animal
    end

    dd = [dd '\' title '.analyzer'];

    ['Saving analyzer file at location:  ' dd]

    save(dd ,'Analyzer')
    
end

