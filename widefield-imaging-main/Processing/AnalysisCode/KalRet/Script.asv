%Script for analyzing kalatsky retinotopy with 2photon data

f1 = f1meanimage;  %Build F1 images (takes the longest)
L = fspecial('gaussian',15,1);  %make spatial filter
bw = ones(size(f1{1}));
[kmap_hor kmap_vert] = processkret(f1,bw,L);  %Make maps to plot, delete L if no smoothing

%%%%%%%%%%%%%%%

figure
imagesc(kmap_hor,[-180 180])
title('Horizontal Retinotopy')
colorbar('SouthOutside')
colormap hsv
truesize

figure
imagesc(kmap_vert,[-180 180])
title('Vertical Retinotopy')
colorbar
colormap hsv
truesize

%%%%%%%%%%%%%%%%%

figure
contour(kmap_hor.*bw,-180:20:180,'r')
hold on
contour(kmap_vert.*bw,-180:20:180,'b')
title('Contour')
axis ij     %'contour' plots inverted

%%%%%%%%%%%%%%%%%%

figure
imagesc(sh)
colormap gray
title('ROI Coverage of Stimulus Area')