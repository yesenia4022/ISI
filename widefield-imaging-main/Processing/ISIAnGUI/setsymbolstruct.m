function setsymbolstruct

%Put all the symbol information into global structure

global symbolInfo Analyzer G_handles

symbolInfo = struct;


Nsym = length(Analyzer.loops.conds{1}.symbol);

for i = 1:Nsym
    
    switch i
        
        case 1
            
            Fsymbol = get(G_handles.primSymbol,'string'); %primary parameter symbol in looper to analyze
            symbolInfo.ID(1) = get(G_handles.primSymbol,'value');  %The index with respect to the looper
            symbolInfo.str{1} = Fsymbol{symbolInfo.ID(1)};  %Selected string
            symbolInfo.domType = get(G_handles.domType,'value');  %Type of domain for primary symbol... .e.g. circular 'Axis'
            
            
        case 2
            
            Fsymbol = get(G_handles.secSymbol,'string'); %secondary symbol
            symbolInfo.ID(2) = get(G_handles.secSymbol,'value');
            symbolInfo.str{2} = Fsymbol{symbolInfo.ID(2)};
            symbolInfo.Collapse(1) = get(G_handles.secCollapse,'value');
            
        case 3
            
            Fsymbol = get(G_handles.tertSymbol,'string'); %tertiary symbol
            symbolInfo.ID(3) = get(G_handles.tertSymbol,'value');
            symbolInfo.str{3} = Fsymbol{symbolInfo.ID(3)};
            symbolInfo.Collapse(2) = get(G_handles.tertCollapse,'value');
            
        case 4
            
            Fsymbol = get(G_handles.quatSymbol,'string'); %tertiary symbol
            symbolInfo.ID(4) = get(G_handles.quatSymbol,'value');
            symbolInfo.str{4} = Fsymbol{symbolInfo.ID(4)};
            symbolInfo.Collapse(3) = get(G_handles.quatCollapse,'value');
            
        case 5
            
            Fsymbol = get(G_handles.quintSymbol,'string'); %tertiary symbol
            symbolInfo.ID(5) = get(G_handles.quintSymbol,'value');
            symbolInfo.str{5} = Fsymbol{symbolInfo.ID(5)};
            symbolInfo.Collapse(4) = get(G_handles.quintCollapse,'value');
            
    end
    
    
end
