function [xpos ypos xsize ysize] = getPosSize

global Analyzer

%Used for coarse retinotopy analysis
screenXcm = Analyzer.M.screenXcm;
screenYcm = Analyzer.M.screenYcm; %these are dependencies for 'formula' below
screenDist = Analyzer.M.screenDist;
x_zoom = getparam('x_zoom');
y_zoom = getparam('y_zoom');
loopform = Analyzer.L.formula;

for j = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp('a',Analyzer.loops.conds{1}.symbol{j})
        aid = j;
    elseif strcmp('b',Analyzer.loops.conds{1}.symbol{j})
        bid = j;
    end
end

nc = getnoconditions;
bflag = stimblank(nc);
if bflag
    nloop = nc-1;
else
    nloop = nc;
end

for c = 1:nloop

    a = Analyzer.loops.conds{c}.val{aid};
    b = Analyzer.loops.conds{c}.val{bid};

    eval([loopform ';']);

    if b == 0
        xpos(c) = x_pos;
        xsize(c) = xP*x_zoom; %units of pixels
        ypos(c) = NaN;
        ysize(c) = NaN;
    else
        ypos(c) = y_pos;
        ysize(c) = yP*y_zoom;
        xpos(c) = NaN;
        xsize(c) = NaN;
    end

end


