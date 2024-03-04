function setGUIlabels

global Analyzer G_handles

Gsetdirectories

conds = getnoconditions;
reps = getnorepeats(1);
set(G_handles.nocond,'string',num2str(conds))
set(G_handles.norep,'string',num2str(reps))
set(G_handles.dirstatus,'string','Loaded')

set(G_handles.setROI,'enable','on')
set(G_handles.process,'enable','on')

set(G_handles.predelay,'string',['predelay=' num2str(getParamVal('predelay',0))])
set(G_handles.postdelay,'string',['postdelay=' num2str(getParamVal('postdelay',0))])
set(G_handles.trialtime,'string',['trialtime=' num2str(getParamVal('stim_time',0))])
    
Nsym = length(Analyzer.loops.conds{1}.symbol);


set(G_handles.primSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
set(G_handles.primSymbol,'value',1)

if Nsym > 1
    set(G_handles.secSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
    set(G_handles.secSymbol,'value',2)
    set(G_handles.secSymbol,'enable','on')
    set(G_handles.secCollapse,'enable','on')
    dom = getdomain(Analyzer.loops.conds{1}.symbol{2});
    domstr{1} = 'mean';
    for i = 2:length(dom)+1
        domstr{i} = dom(i-1);
    end
    set(G_handles.secCollapse,'string',domstr)
else
    set(G_handles.secSymbol,'enable','off')
    set(G_handles.tertSymbol,'enable','off')
    set(G_handles.secCollapse,'enable','off')
    set(G_handles.tertCollapse,'enable','off')
end

if Nsym > 2
    set(G_handles.tertSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
    set(G_handles.tertSymbol,'value',3)
    set(G_handles.tertSymbol,'enable','on')
    set(G_handles.tertCollapse,'enable','on')
    dom = getdomain(Analyzer.loops.conds{1}.symbol{3});
    for i = 1:length(dom)
        domstr{i} = dom(i);
    end
    set(G_handles.tertCollapse,'string',domstr)
else
    set(G_handles.tertSymbol,'enable','off')
    set(G_handles.tertCollapse,'enable','off')
end


% for i = 1:Nsym
%     
%     switch i
%         
%         case 1
%             
%             set(G_handles.primSymbol,'string',Analyzer.loops.conds{1}.symbol);  %Set the strings in drop down menu
%             set(G_handles.primSymbol,'value',i)            
%             
%         case 2
%             
%             set(G_handles.secSymbol,'string',Analyzer.loops.conds{1}.symbol);
%             set(G_handles.secSymbol,'value',i)
%             set(G_handles.secSymbol,'enable','on')
%             set(G_handles.secCollapse,'enable','on')
%             
%         case 3
%             
%             set(G_handles.tertSymbol,'string',Analyzer.loops.conds{1}.symbol);
%             set(G_handles.tertSymbol,'value',i)
%             set(G_handles.tertSymbol,'enable','on')
%             set(G_handles.tertCollapse,'enable','on')
%             
%         case 4
%             
%             set(G_handles.quatSymbol,'string',Analyzer.loops.conds{1}.symbol);
%             set(G_handles.quatSymbol,'value',i)
%             set(G_handles.quatSymbol,'enable','on')
%             set(G_handles.quatCollapse,'enable','on')
%             
%         case 5
%             
%             set(G_handles.quintSymbol,'string',Analyzer.loops.conds{1}.symbol);
%             set(G_handles.quintSymbol,'value',i)
%             set(G_handles.quintSymbol,'enable','on')
%             set(G_handles.quintCollapse,'enable','on')            
%             
%     end
%    
%     
% end
% 
% %Turn off the unused fields
for i = Nsym+1:5
    
    switch i
        
        case 2
            
            set(G_handles.secSymbol,'enable','off')
            set(G_handles.secCollapse,'enable','off')            
            
        case 3
            
            set(G_handles.tertSymbol,'enable','off')
            set(G_handles.tertCollapse,'enable','off')
            
        case 4
            
            set(G_handles.quatSymbol,'enable','off')
            set(G_handles.quatCollapse,'enable','off')
            
        case 5
            
            set(G_handles.quintSymbol,'enable','off')
            set(G_handles.quintCollapse,'enable','off')
            
    end
end

