function [mat xdom ydom] = smoothscatter(x,y,sigx,sigy,varargin)

%Make grid

N = 200;

if ~isempty(varargin)
    xlimits = varargin{1};
    ylimits = varargin{2};
else     
    if min(y) > 0
        ylimits = [0 prctile(y,99.9)];
    else
        ylimits = [prctile(y,.1) prctile(y,99.9)];
    end
    
    if min(x) > 0
        xlimits = [0 prctile(x,99.9)];
    else
        xlimits = [prctile(x,.1) prctile(x,99.9)];
    end
    
end

xdom = linspace(xlimits(1),xlimits(2),N);
ydom = linspace(ylimits(1),ylimits(2),N);

dx = xdom(2)-xdom(1);
dy = ydom(2)-ydom(1);

[dum id] = min(abs(xdom-0));
xdom = xdom-xdom(id);
[dum id] = min(abs(ydom-0));
ydom = ydom-ydom(id); %make sure point is centered on zero

mat = zeros(length(ydom),length(xdom));

for i = 1:length(ydom);
    id = find(y>ydom(i)-dy/2 & y<ydom(i)+dy/2);
    h = hist(x(id),xdom);
    mat(i,:) = h;
end

[xdomMat ydomMat] = meshgrid(xdom-mean(xdom),ydom-mean(ydom));
G = exp(-xdomMat.^2/(2*sigx^2)) .* exp(-ydomMat.^2/(2*sigy^2));

%mat = ifft2(abs(fft2(G)).*fft2(mat));
mat = conv2(mat,G,'same');