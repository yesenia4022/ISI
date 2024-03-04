function bwV = getVesselMask

global bw

im = getTrialMean([0 4000],1,1);
%%
R = 20;
h = ones(R,R)/(R^2);
mu = filter2(h,im);  %Compute local mean
imout = im-mu;

st = sqrt(filter2(h,imout.^2));  %local stddev

imout = imout./st;  %Compute local Z value

id = find(isnan(imout));
imout(id) = 0;

id = find(~bw);
imout(id) = 0;
%%
thresh = .8;
imout2 = (sign(imout-thresh)+1)/2;

mmorph = 1;
SE = strel('disk',mmorph,0);
imout2 = imopen(imout2,SE);

%%
minsize = 50;
%Get rid of cells that are smaller than minsize
if minsize > 0
    
    celllabel = bwlabel(imout2);
    cellid = unique(celllabel);
    for j = 2:length(cellid)
        id = find(cellid(j) == celllabel);
        length(id)
        if length(id) < minsize
            imout2(id) = 0;
        end
    end
    
end
    
%%
bwV = imout2;
figure,imagesc(imout2), colormap gray