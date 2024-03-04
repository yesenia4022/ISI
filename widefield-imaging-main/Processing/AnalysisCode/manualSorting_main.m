ManualCellSorter(im)

[bwCell1 bwCell2] = MakeCellMask(14,.8,5,9);

dim = size(bwCell2);
posallvec = (posall(:,1)-1)*dim(1) + posall(:,2);

celllabel = bwlabel(bwCell2);
cellid = unique(celllabel);
for i = 2:length(cellid)
    id = find(cellid(i) == celllabel);
    c = intersect(id,posallvec);
    if isempty(c)
        bwCell2(id) = 0;
    end
end

