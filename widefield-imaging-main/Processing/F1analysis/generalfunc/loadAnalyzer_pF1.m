function loadAnalyzer(ue)

global anadir Analyzer
anadir
Anim = anadir(end-3:end-1);
fname = [Anim '_' ue '.analyzer'];
path = [anadir fname]

load(path,'-mat')

%Get the syncInfo and put it in the Analyzer structure for convenience

i = 1;
while exist(['syncInfo' num2str(i)])    
   eval(['syncInfoAll{i} = syncInfo' num2str(i) ';']) 
   i = i+1;   
end

Analyzer.syncInfo = syncInfoAll;