function plotF0(f0,bw,hh)

N = sqrt(length(f0));


for i = 1:length(f0)
    
    if ~isempty(hh)
        f0{i} = ifft2(fft2(f0{i}).*abs(fft2(hh)));
    end
    
    %f0dum(find(~bw)) = mean(f0dum(:));
    
    midum = prctile(f0{i}(find(bw)),1);
    madum = prctile(f0{i}(find(bw)),99);
    
    if i == 1
        mi = midum;
        ma = madum;
    end
    
    ma = max([madum ma]);
    mi = min([midum mi]);
    
end


[y x] = find(bw);
ywin = min(y):max(y);
xwin = min(x):max(x);

figure
for i = 1:length(f0)

    %subplot(ceil(N),floor(N),i)
    subplot(round(sqrt(length(f0))),ceil(sqrt(length(f0))),i)
    dum = f0{i};
    dum(~bw) = mi;
    imagesc(dum(ywin,xwin),[mi ma]), axis image
    %set(gca,'Xtick',[],'Ytick',[]), 
    colorbar
    
    try
        hold on,
        plot([112 112],[68 68],'xr')
    catch
    end
    

end

colormap gray