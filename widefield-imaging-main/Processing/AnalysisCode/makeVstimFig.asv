N = 256;

[x y] = meshgrid(0:N-1,0:N-1);

mid = ceil(N/2);
r = sqrt((x-mid).^2 + (y-mid).^2);
id = find(r<=mid);
Mask = zeros(N,N);
Mask(id) = 1;

Ncycles = [3 8 4];

ori = [0 30 100]*pi/180;
gain{1} = [1 1 1]; gain{2} = [1 1 1]; gain{3} = [1 1 1];
%gain{1} = [1 -.18 -.14]; gain{2} = [-1 .52 -.027]; gain{3} = [.21 -.27 1]; %LMS cone

figure
for i = 1:3

    x2 = (x/N)*Ncycles(i)*2*pi;
    y2 = (y/N)*Ncycles(i)*2*pi;
    xp = x2*cos(ori(i)) + y2*sin(ori(i));

    im(:,:,1) = (sin(xp).*Mask*gain{i}(1)+1)/2;
    im(:,:,2) = (sin(xp).*Mask*gain{i}(2)+1)/2;
    im(:,:,3) = (sin(xp).*Mask*gain{i}(3)+1)/2;
    
    subplot(1,3,i)
    image(im), colormap gray
    
    axis image
    axis off
end

