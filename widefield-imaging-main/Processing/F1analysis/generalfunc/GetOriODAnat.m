function GetOriODAnat(fw,orirng,odrng)

%Need to plot in processF0 and have anatomy image loaded from overlay GUI
%fw is intensity of the colormap relative to the anatomy [0 1]

global f0m1 f0m2 imstate bw

reset_imstate_F0;

hh = makeFilt;

%Get orientation map
funcmap1 = GprocessOri(f0m1,hh);  %Channel 1
funcmap2 = GprocessOri(f0m2,hh);
funcmap = (funcmap1+funcmap2)/2;
orimap = angle(funcmap);
orimap = (orimap+pi*(1-sign(orimap)))/2*180/pi;

%Get ocular dominance map
funcmap1 = GprocessOcdom(f0m1,hh);
funcmap2 = GprocessOcdom(f0m2,hh);
funcmap = (funcmap1+funcmap2)/2;
funcmap(find(funcmap<prctile(funcmap(:),.5))) = prctile(funcmap(:),.5);
funcmap(find(funcmap>prctile(funcmap(:),99.5))) = prctile(funcmap(:),99.5);
odmap = funcmap;

odmap(find(1-bw)) = NaN;
orimap(find(1-bw)) = NaN;

odmap = fliplr(flipud(odmap'));
orimap = fliplr(flipud(orimap'));



%%

try
    imanat = imstate.imanat;
catch
    'load anatomy'
    adfsn 
end
    
%imanat(find(1-bw)) = NaN;

imanat(1,1) = 0;  %in case they are all ones, the next operations don't make sense

imanat(find(imanat<prctile(imanat(:),.5))) = prctile(imanat(:),.5);
imanat(find(imanat>prctile(imanat(:),99.5))) = prctile(imanat(:),99.5);

imanat = imanat-min(imanat(:));
imanat = imanat/max(imanat(:));
imanat = round(imanat*63+1);

id = find(isnan(imanat));
imanat(id) = nanmedian(imanat(:));

imanat = fliplr(flipud(imanat'));

%% Make ORI image

imfunc = orimap;
imfunc = imfunc-min(imfunc(:));
imfunc = imfunc/max(imfunc(:));
imfunc = round(imfunc*63+1);

id = find(isnan(imfunc));
imfunc(id) = nanmedian(imfunc(:));

figure, grayid = gray; close;

figure, hsvid = hsv; close
dim = size(imfunc);

aw = 1-fw;
for i = 1:dim(1)
    for j = 1:dim(2)
        imout(i,j,:) = fw*hsvid(imfunc(i,j),:) + aw*grayid(imanat(i,j),:);       
    end
end
imout = imout/max(imout(:));

figure
subplot(1,2,1)
image(imout)
colorbar
colormap hsv
colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})
title('orientation map')
axis image

subplot(1,2,2)
image(imout)
colorbar
colormap hsv
colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})
title('orientation map; od contours')
axis image

hold on, contour(odmap,[0 0],'k','LineWidth',1)

%% Make OD image

imfunc = odmap;
imfunc = imfunc-min(imfunc(:));
imfunc = imfunc/max(imfunc(:));
imfunc = round(imfunc*63+1);

id = find(isnan(imfunc));
imfunc(id) = nanmedian(imfunc(:));

figure, jetid = jet; close 
dim = size(imfunc);

aw = 1-fw;
for i = 1:dim(1)
    for j = 1:dim(2)
        imout(i,j,:) = fw*jetid(imfunc(i,j),:) + aw*grayid(imanat(i,j),:);       
    end
end
imout = imout/max(imout(:));
 

figure
subplot(1,2,1)
image(imout)
colormap jet
dum = linspace(min(odmap(:)),max(odmap(:)),5);
for i = 1:length(dum)
    cbardom{i} = dum(i);
end
colorbar('YTick',[1 16:16:64],'YTickLabel',cbardom)
title('od map (R-L)')
axis image

subplot(1,2,2)
image(imout)
colormap jet
dum = linspace(min(odmap(:)),max(odmap(:)),5);
for i = 1:length(dum)
    cbardom{i} = dum(i);
end
colorbar('YTick',[1 16:16:64],'YTickLabel',cbardom)
title('od map; ori contours')
axis image

hold on
contour(orimap/180,[0:45:180]/180,'k')

%% Anatomy with intersection points

id = find(orimap>orirng(1) & orimap<orirng(2) & odmap>odrng(1) & odmap<odrng(2));
bounds = zeros(size(orimap));
bounds(id) = 1;

figure, 

subplot(1,2,1)
imagesc(imanat), colormap gray
hold on 
contour(bounds,[.5 .5],'r')
axis image
 
scalemask = bounds*.8 + .2;
subplot(1,2,2)
imagesc(imanat.*scalemask)
hold on 
contour(bounds,[.5 .5],'r')
axis image
colormap gray
axis image
title('Target region highlighted')

function hh = makeFilt

global F0handles f0m1

togstateHP = get(F0handles.HPflag,'Value');
togstateLP = get(F0handles.LPflag,'Value');

sizeIM = size(f0m1{1});
if togstateHP == 1
    Hwidth = str2double(get(F0handles.Hwidth,'string'));
    ind = get(F0handles.HPWind,'value');
    
    switch ind
        case 1
            H = -fspecial('gaussian',sizeIM,Hwidth);
            H(round(sizeIM(1)/2),round(sizeIM(2)/2)) = 1+H(round(sizeIM(1)/2),round(sizeIM(2)/2));
        case 2
            H = hann(Hwidth)*hann(Hwidth)';
            H = -H./sum(H(:));
            H(round(Hwidth/2),round(Hwidth/2)) = 1+H(round(Hwidth/2),round(Hwidth/2));
            Hsize = length(H(1,:));
            marginR = (sizeIM(1)-Hsize)/2;
            marginC = (sizeIM(2)-Hsize)/2;
            H = [zeros(floor(marginR),sizeIM(2)); [zeros(Hsize,floor(marginC)) H zeros(Hsize,ceil(marginC))]; zeros(ceil(marginR),sizeIM(2))];
        case 3
            H = -fspecial('disk',round(Hwidth/2));
            Hsize = length(H(1,:));  %~=Hwidth
            H(round(Hsize/2),round(Hsize/2)) = 1+H(round(Hsize/2),round(Hsize/2));
            marginR = (sizeIM(1)-Hsize)/2;
            marginC = (sizeIM(2)-Hsize)/2;
            H = [zeros(floor(marginR),sizeIM(2)); [zeros(Hsize,floor(marginC)) H zeros(Hsize,ceil(marginC))]; zeros(ceil(marginR),sizeIM(2))];
            
    end
    if togstateLP == 0
        hh = ifft2(abs(fft2(H)));   %Eliminate phase information
    end
end

if togstateLP == 1
    Lwidth = str2double(get(F0handles.Lwidth,'string'));
    ind = get(F0handles.LPWind,'value');
    
    switch ind
        case 1
            L = fspecial('gaussian',sizeIM,Lwidth);
        case 2
            L = hann(Lwidth)*hann(Lwidth)';
            L = L./sum(L(:));
            Lsize = length(L(1,:));
            marginR = (sizeIM(1)-Lsize)/2;
            marginC = (sizeIM(2)-Lsize)/2;
            L = [zeros(floor(marginR),sizeIM(2)); [zeros(Lsize,floor(marginC)) L zeros(Lsize,ceil(marginC))]; zeros(ceil(marginR),sizeIM(2))];
        case 3
            L = fspecial('disk',round(Lwidth/2));
            Lsize = length(L(1,:));
            marginR = (sizeIM(1)-Lsize)/2;
            marginC = (sizeIM(2)-Lsize)/2;
            L = [zeros(floor(marginR),sizeIM(2)); [zeros(Lsize,floor(marginC)) L zeros(Lsize,ceil(marginC))]; zeros(ceil(marginR),sizeIM(2))];
    end
    if togstateHP == 0
        hh = ifft2(abs(fft2(L)));   %Eliminate phase information
    else
        hh = ifft2(abs(fft2(L).*fft2(H)));   %Take mag because phase gives a slight shift.
    end
end


if ~or(togstateLP,togstateHP)
    hh = [];
end
