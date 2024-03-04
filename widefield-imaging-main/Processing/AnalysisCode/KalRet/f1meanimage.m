function f1m = f1meanimage

global f1m
% Compute mean f1 across all conditions and repeats


nc = getnoconditions;
f1 = cell(1,nc);

for c = 1:nc
    f1{c} = f1image(c);
end

% Now average all the repeats

for c = 1:nc
    img = f1{c}{1};
    nr = length(f1{c});
    for(r=2:nr)
        img = img+f1{c}{r};
    end
    img = img/nr;
    f1m{c} = img;  %% Compute mean image
end

