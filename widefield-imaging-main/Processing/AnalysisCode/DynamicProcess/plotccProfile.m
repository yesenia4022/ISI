function plotccProfile(cc,bwCell1)

masklabel = bwlabel(bwCell1);

celldom = unique(masklabel);
celldom = celldom(1:end);
Ncell = length(celldom);

mask = masklabel./masklabel;

figure
for p = 1:Ncell
    
    plotter = zeros(size(bwCell1));
 
    for k = 1:Ncell
        
        idcell = find(masklabel(:) == celldom(k));

        plotter(idcell) = mean(cc(p,k,:));

    end
    
    subplot(ceil(sqrt(Ncell)),ceil(sqrt(Ncell)),p)
    imagesc(plotter,'AlphaData',mask,[-1 1])
    axis off
    
end