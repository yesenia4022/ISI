function [f1m signals] = Gf1meanimage(varargin)

global Analyzer f1 GF1_handles

% Compute mean f1 across all conditions and repeats

nc = length(Analyzer.loops.conds);
bflag = 0;
if strcmp(Analyzer.loops.conds{end}.symbol{1},'blank');
    bflag = 1;
end

f1 = cell(1,nc);
sig1 = cell(1,nc);

set(GF1_handles.status,'string','0%'), drawnow
for c=1:(nc-bflag)
    [f1{c} sig1{c}] = Gf1image(c,varargin);
    set(GF1_handles.status,'string',[num2str(round((c/(nc-bflag))*100)) '%']), drawnow
end

% Now average all the repeats

if length(varargin) == 2
    for c = 1:(nc-bflag)
        nr = length(f1{c});
        sig2 = addtrunc(sig1{c},nr); %sig2 is a matrix where each row is a pixel (mean is already subtracted)
        sig2 = sig2./nr;
        signals{c} = sig2;
    end
else
    signals = 0;
end


for c=1:(nc-bflag)
    img = f1{c}{1};
    nr = length(f1{c});
    for r=2:nr
        img = img+f1{c}{r};
    end
    img = img/nr;
    f1m{c} = img;  %% Compute mean image
end

set(GF1_handles.status,'string','Done'), drawnow


function y = addtrunc(x,nr)
%Truncates all signals in x to the length of shortest one, and then adds them.
for i = 1:nr
    N(i) = length(x{i}(1,:));
end
shortest = min(N);  %Length of shortest repeat

y = zeros(length(x{1}(:,1)),shortest); %(No. of pixels) x (length of shortest repeat)
for i = 1:nr
    y = y + x{i}(:,1:shortest);
end
 
