function loadAnalyzer(ue)

global anadir Analyzer G_handles
anadir
Anim = get(G_handles.loadana,'string');
fname = [Anim '_' ue];
path = [anadir fname '.analyzer']

load(path,'-mat')

%Get the syncInfo and put it in the Analyzer structure for convenience

i = 1;
while exist(['syncInfo' num2str(i)])    
   eval(['syncInfoAll{i} = syncInfo' num2str(i) ';']) 
   i = i+1;   
end

Analyzer.syncInfo = syncInfoAll;