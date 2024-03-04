function [f0 ori cont] = GetOrivsContrast(f0,hh)

global Analyzer
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping params
for idsym = 1:Nsym
    if strcmp('ori',Analyzer.loops.conds{1}.symbol{idsym})
        oriid = idsym;
    elseif strcmp('contrast',Analyzer.loops.conds{1}.symbol{idsym})
        contid = idsym;
    end
end

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    f0blank = f0{end};
    f0(end) = [];
end

for i = 1:length(f0)
    ori(i) = Analyzer.loops.conds{i}.val{oriid};
    cont(i) = Analyzer.loops.conds{i}.val{contid};
end


%if a filter exists, use it...
if ~isempty(hh)
    id = find(isnan(orimap));
    orimap(id) = 0;
    orimap = ifft2(abs(fft2(hh)).*fft2(orimap));    
end

