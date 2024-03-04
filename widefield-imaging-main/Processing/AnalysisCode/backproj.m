function im = backproj(kern,rotdom)

%1st dimension of kern should be position, and 2nd should be rotation

dim = size(kern); 
%barW = 5;  %bar width
%sig = barW/3;
%imW = barW*dim(1);  %image width
sig = .5;
imW = length(kern(:,1));

dom = -floor(imW/2):ceil(imW/2)-1;
[x y] = meshgrid(dom,fliplr(dom));

im = zeros(imW,imW);

for i = 1:dim(2)
    xp = x*cos(rotdom(i)*pi/180) + y*sin(rotdom(i)*pi/180);
    yp = -x*sin(rotdom(i)*pi/180) + y*cos(rotdom(i)*pi/180);

    
    for j = 1:dim(1)
        
%         Llim = (j-1)*barW+1;
%         Rlim = j*barW;
%         
%         id = find(xp>=dom(Llim) & xp<=dom(Rlim));
% 
%         imbar = zeros(imW,imW);
%         imbar(id) = kern(i,j);
        
        %xc = dom(round(barW*j - barW/2));
        yc = dom(j);
        imbar = kern(j,i)*exp(-(xp-yc).^2/(2*sig^2));

        im = im + imbar;
        
        %imagesc(im)
        %hold on
        %pause(.4)

    end
end

r = sqrt(x.^2 + y.^2);
id = find(r>dom(end));
im(id) = NaN;

        