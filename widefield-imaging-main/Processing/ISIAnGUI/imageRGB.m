function imageRGB(im,idx)
%idx is 1,2,or 3, corresponding to R,G,or B
%plot R G B channel

dim = size(im);
im = (im-min(im(:)));
im = im/max(im(:));
plotter = zeros(dim(1),dim(2),3);
plotter(:,:,idx) = im; 

image(plotter)       %Load and plot frame