function plotF0dum(f0,bw,hh)

N = sqrt(length(f0));


load('C:\Users\ISI\Desktop\OG figs\imstate_myb_001_000')
%load('C:\Users\ISI\Desktop\OG figs\na2\imstate_na2_001_001')
%bw = imstate.areaBounds;

for i = 1:length(f0)
    if ~isempty(hh)
        f0{i} = ifft2(fft2(f0{i}).*abs(fft2(hh)));
    end 
end


for i = 1:length(f0)

    %f0norm{i} = f0{i}-f0{end};
    f0norm{i} = f0{i};
    
    midum = prctile(f0norm{i}(find(bw)),50);
    madum = prctile(f0norm{i}(find(bw)),80);
    
    if i == 1
        mi = midum;
        ma = madum;
    end
    
    ma = max([madum ma]);
    mi = min([midum mi]);
    
end

%mi = 0;

%%
[y x] = find(bw);
ywin = min(y):max(y);
xwin = min(x):max(x);

%refpt = [157 113]; %[y x]

refpt = [319 287] - [ywin(1) xwin(1)];

refpt = 225 - refpt;

refpt = [125 100]

mapOverlayBit = 1

if mapOverlayBit
    load('C:\Users\ISI\Desktop\OG figs\imstate_myb_001_000')
    %load('C:\Users\ISI\Desktop\OG figs\na2\imstate_na2_001_001')
    
    hor_map = imstate.fmaps{1}(ywin,xwin);
    vert_map = imstate.fmaps{2}(ywin,xwin);
    Vareas = round(imstate.areaBounds(ywin,xwin));
    sigMag = imstate.sigMag(ywin,xwin);
    %sigMag = log(sigMag);

    l = prctile(sigMag(:),5); 
    h = prctile(sigMag(:),90);
    sigMag(find(sigMag<l)) = l;
    sigMag(find(sigMag>h)) = h;
    
    sigMag = sigMag-min(sigMag(:))+.0001;
    sigMag = sigMag/max(sigMag(:));
    
    
    hor_map = imrotate(hor_map,180);
    vert_map = imrotate(vert_map,180);
    Vareas = imrotate(Vareas,180);
    sigMag = imrotate(sigMag,180);
    
    figure,
    subplot(1,2,1)
    imagesc(hor_map,'AlphaData',sigMag,[-50 50]), colormap jet
    hold on ,
    contour(Vareas,[.5 .5], 'k','LineWidth',2);
    hold on,
    plot(refpt(2),refpt(1),'xk','MarkerSize',10,'LineWidth',2)
    axis image
    colorbar
    
    subplot(1,2,2)
    imagesc(vert_map,'AlphaData',sigMag,[-50 50]), colormap jet
    hold on ,
    contour(Vareas,[.5 .5], 'k','LineWidth',2);
    hold on,
    plot(refpt(2),refpt(1),'xk','MarkerSize',10,'LineWidth',2)
    axis image
    colorbar
end


figure
for i = 1:length(f0)

    dum = f0norm{i};
    dum(~bw) = mi;
    f0map = dum(ywin,xwin);
    f0map = imrotate(f0map,180);
    
    subplot(round(sqrt(length(f0))),ceil(sqrt(length(f0))),i)
    imagesc(f0map,[mi ma]), axis image
    %set(gca,'Xtick',[],'Ytick',[]),
    colorbar    
    
    if mapOverlayBit
        hold on ,
        contour(Vareas,[.5 .5], 'y','LineWidth',2);
        hold on,
        plot(refpt(2),refpt(1),'xr','MarkerSize',10,'LineWidth',2)
        axis image
        colorbar
    end    
    
    try
        hold on,
        plot(refpt(2),refpt(1),'xr')
    catch
    end

end
colormap gray

%%


%%
