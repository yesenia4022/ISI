
load('E:\ISIdata\se6\grabs\se6_U1_vasculature___grab_se6_001_000_24_Aug_2018_08_26_39.mat')
load('E:\ISIdata\se7\grabs\se7_vasculature_after_U0____grab_se7_000_000_24_Aug_2018_11_37_56.mat')
figure
%vasc = grab.img(1:606,1:808);
vasc = grab.img-1000;
imagesc(vasc)
%colormap gray
hold on

vert = kmap_vert-min(kmap_vert(:));
vert = vert.*28;
h = imagesc(vert);
%colormap jet
I = .3*ones(size(vasc));
set(h, 'AlphaData', I)

hor = kmap_hor-min(kmap_hor(:));
hor = hor.*35;
h = imagesc(hor);
%colormap jet
I = .5*ones(size(vasc));
set(h, 'AlphaData', I)