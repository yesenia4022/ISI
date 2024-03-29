function updateLstate

%Ian Nauhaus

global Lstate GUIhandles

Lstate.reps = str2num(get(GUIhandles.looper.repeats,'string'));
Lstate.blankPeriod = str2num(get(GUIhandles.looper.blankPeriod,'string'));
Lstate.rand = get(GUIhandles.looper.randomflag,'value');
Lstate.blank = get(GUIhandles.looper.blankflag,'value');

Ldum{1} = {[get(GUIhandles.looper.symbol1,'string')] [get(GUIhandles.looper.valvec1,'string')]};
Ldum{2} = {[get(GUIhandles.looper.symbol2,'string')] [get(GUIhandles.looper.valvec2,'string')]};
Ldum{3} = {[get(GUIhandles.looper.symbol3,'string')] [get(GUIhandles.looper.valvec3,'string')]};
Ldum{4} = {[get(GUIhandles.looper.symbol4,'string')] [get(GUIhandles.looper.valvec4,'string')]};
Ldum{5} = {[get(GUIhandles.looper.symbol5,'string')] [get(GUIhandles.looper.valvec5,'string')]};


%Get rid of blank rows...
Lstate.param = cell(1,1);  %initialize
k = 1;
for i = 1:length(Ldum)
    if ~isempty(Ldum{i}{1})
        Lstate.param{k} = Ldum{i};        
        k = k+1;
    end
end

Lstate.formula = get(GUIhandles.looper.formula,'string');