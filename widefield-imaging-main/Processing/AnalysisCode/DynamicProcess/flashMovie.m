function flashMovie(kern,bwCell1)



masklabel = bwlabel(bwCell1);

celldom = unique(masklabel);
celldom = celldom(1:end);
Ncell = length(celldom);

mask = masklabel./masklabel;

%make movie

plotter = zeros(length(bwCell1(:,1,1)),length(bwCell1(1,:,1)),length(kern{1}));
for t = 1:length(kern{1})
    plotterdum = zeros(size(bwCell1));
    for p = 2:Ncell

        idcell = find(masklabel(:) == celldom(p));

        plotterdum(idcell) = kern{p}(t);

    end
    plotter(:,:,t) = plotterdum;
end

figure
for t = 1:length(plotter(1,1,:))
    
    imagesc(plotter(:,:,t),[min(plotter(:)) max(plotter(:))]), colormap gray
    pause(.2)
    
end
