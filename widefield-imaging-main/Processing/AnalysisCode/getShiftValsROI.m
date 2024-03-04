function [mbest nbest] = getShiftValsROI(im,temp,bwx,bwy)

global processF0_handles

dim = size(im);

W = str2double(get(processF0_handles.searchRange,'String'));
 
bwid = dim(1)*(bwx-1) + bwy;

shiftdom = -W:W;

i = 0; 
CC = zeros(length(shiftdom),length(shiftdom));
for dy = -W:W
    i = i+1;
    j = 0;
    for dx = -W:W
        j = j+1;
        bwx2 = bwx+dx;
        bwy2 = bwy+dy;
        
        idIB = find(bwx2<=dim(1) & bwx2>0 & bwy2<=dim(2) & bwy2>0);  %In bounds
        
        bwid2 = dim(1)*(bwx2-1) + bwy2;
        
        Rdum = corrcoef(temp(bwid(idIB)),im(bwid2(idIB)));
        
        CC(i,j) = Rdum(1,2);
        
    end
end

[mbest nbest] = find(CC == max(CC(:)));
mbest = shiftdom(mbest(1));
nbest = shiftdom(nbest(1));


