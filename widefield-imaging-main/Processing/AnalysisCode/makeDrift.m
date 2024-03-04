function stack = makeDrift(ori,Ncycles)

ori = ori*pi/180;

N = 256;

[x y] = meshgrid(0:N-1,0:N-1);

mid = ceil(N/2);
r = sqrt((x-mid).^2 + (y-mid).^2);
id = find(r<=mid);
Mask = zeros(N,N);
Mask(id) = 1;

x = (x/N)*Ncycles*2*pi;
y = (y/N)*Ncycles*2*pi;

xp = x*cos(ori) + y*sin(ori);

res = 15;

figure(1)
for i = 1:res*Ncycles
   phase = (i-1)*2*pi/res; 
   im = cos(xp+phase).*Mask;
   imagesc(im), colormap gray;
   set(gca,'Xtick',[],'Ytick',[])
   axis image;
   stack(i) = getframe;
end

movie2avi(stack,'C:\Documents and Settings\SNLC\Desktop\grating','fps',25)