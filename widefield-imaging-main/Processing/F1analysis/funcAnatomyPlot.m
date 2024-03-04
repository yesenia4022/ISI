function funcAnatomyPlot

global imstate

aw = 1-imstate.intRatio;  %anatomy weight of image (scalar)
fw = imstate.intRatio;  %anatomy weight of image (scalar)
imanat = imstate.imanat;
imfunc = imstate.imfunc;

if unique(imanat(:)) ~= 1
    imanat(find(imanat<prctile(imanat(:),.5))) = prctile(imanat(:),.5);
    imanat(find(imanat>prctile(imanat(:),99.5))) = prctile(imanat(:),99.5);
    imanat = imanat-min(imanat(:));
    imanat = imanat/max(imanat(:));
    imanat = round(imanat*63+1);
end



imfunc = imfunc-min(imfunc(:));
imfunc = imfunc/max(imfunc(:));
imfunc = round(imfunc*63+1);


grayid = gray;
hsvid = hsv;
%%
mag = imstate.mag.*imstate.bw;
%mag = log(mag);
% mag = medfilt2(mag,[3 3]);

% h = fspecial('gaussian',size(mag),.5);
% h = abs(fft2(h));
% magf = ifft2(h.*fft2(mag));

%mag = magf.^1;
if unique(imanat(:)) ~= 1
    mag = mag-min(mag(:));
    mag = mag/max(mag(:));
end

% 
% thresh = .12;
% id = find(mag(:)<thresh);
% mag(id) = 0;
% 
% 
% figure(29),imagesc(imfunc,'AlphaData',mag), colormap hsv
%%


dim = size(imfunc);

for i = 1:dim(1)
    for j = 1:dim(2)
        if isnan(imfunc(i,j))
            imfunc(i,j) = 1;
        end    
        imout(i,j,:) = fw*mag(i,j)*hsvid(imfunc(i,j),:) + aw*grayid(imanat(i,j),:);       
    end
end
imout = imout/max(imout(:));
 

figure(86)
image(imout)

axis image

