function [xpos ypos omap BWmap] = plotRandPosMap_xy2(blockSize,va_thresh)

%Ian Nauhaus

global kernelsIm G_RChandles ACQinfo

%Get the time domain
eval(['kernDel = ' get(G_RChandles.kernelLength,'string')  ';']);
tauL = kernDel(2)-kernDel(1); %ms
acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine; 
Ntau = round(tauL/acqPeriod)+1;
taudom = (0:Ntau-1)*acqPeriod + kernDel(1);  

%Compute the orimap
dim = size(kernelsIm);
dim2 = size(kernelsIm{1,1,1}); 

t1 = 100; t2 = 600;
[dum id1] = min(abs(taudom-t1));
[dum id2] = min(abs(taudom-t2));

xdom = 0:(dim(2)-1);
xdom = xdom-mean(xdom);

dori = 180/dim(1);
oridom = (0:dori:(180-dori)) + getparam('ori');


% DS = 0;  %downsample
% if DS
%     for ori = 1:dim(1)
%         Tens{ori} = (Tens{ori}(:,:,1:2:end) + Tens{ori}(:,:,2:2:end))/2;
%     end
%     xdom = xdom(2:2:end); xdom = xdom-mean(xdom);
%     dim(2) = length(xdom);
%     
%     for i = 1:dim(1)/2
%         Tens2{i} = Tens{i*2-1} + Tens{i*2};
%     end
%     Tens = Tens2;  clear Tens2
%     oridom = (oridom(1:2:end) + oridom(2:2:end))/2;
%     dim(1) = length(oridom);
%     
% end



h = fspecial('gaussian', size(kernelsIm{1,1,1,1}(:,:,1)), 1);
h = abs(fft2(h));

%Make a spatial tuning curve for each orientation
clear v
Nblocks = ACQinfo.linesPerFrame/blockSize;
for ori = 1:dim(1)
    for bw = 1:2
        for x = 1:dim(2)            
            
            imdum = mean(kernelsIm{ori,x,1,bw}(:,:,id1:id2),3);
            
            %imdum = ifft2(fft2(imdum).*h);
            
            for q = 1:Nblocks  %get position/orientation kernel for each block
                for p = 1:Nblocks
                    dum = imdum(1+(q-1)*blockSize:q*blockSize,1+(p-1)*blockSize:p*blockSize);
                    v{q,p}(bw,x,ori) = mean(dum(:));
                end
            end            
            
        end
    end
end

clear BW
figure
k = 1;
for q = 1:Nblocks  %get position/orientation kernel for each block
    for p = 1:Nblocks
        dum1 = v{q,p}(1,:,:); dum2 = v{q,p}(2,:,:);
        BWmap(q,p) = log(var(dum1(:))./var(dum2(:)));
        v{q,p} = squeeze(mean(v{q,p},1));
        
        subplot(Nblocks,Nblocks,k)
        imagesc(v{q,p})
        k = k+1;
    end
end

figure,imagesc(BWmap), axis image
%%

degperpix = getparam('x_size')/(length(v{1,1}(:,1))-1);  %can't use xdom because of the padding

clear RF xpos ypos omap
k = 1;
%figure
for q = 1:Nblocks
    q
    for p = 1:Nblocks
        
        vdum = v{q,p};
        
        thresh1 = prctile(vdum(:),2);
        vdum = phi(vdum-thresh1);
        
        RF{q,p} = iradon(vdum,oridom,'none','linear',1,length(v{q,p}(:,1)));
        
        
        
        %RF{q,p}(:,end/2:end) = median(RF{q,p}(:));
        %RF{q,p}(:,1:end/2) = median(RF{q,p}(:));
        
%         xdom = (0:length(RF(1,:))-1)*degperpix; %degrees
%         ydom = (0:length(RF(:,1))-1)*degperpix;   %don't use 'y_size'

        thresh2 = prctile(RF{q,p}(:),20) + (max(RF{q,p}(:))-prctile(RF{q,p}(:),20))*.5;

        RF{q,p} = phi(RF{q,p}-thresh2);

        RF{q,p} = RF{q,p}/max(RF{q,p}(:));
        
%         figure(98)
%         hold on, contour(RF{q,p},[.5 .5], 'k')
        
        otc = max(vdum);

        [dum omap(q,p)] = orifind(otc',oridom);

        [param ffit varaccount] = Gaussfit2Drot(double(RF{q,p}),double(omap(q,p)));

        if varaccount > va_thresh
            xpos(q,p) = param(2)*degperpix;
            ypos(q,p) = param(1)*degperpix;
        else
            xpos(q,p) = NaN;
            ypos(q,p) = NaN;
        end

%         subplot(Nblocks,Nblocks,k)
%         imagesc(RF{q,p}), drawnow
%         axis square
%         axis off
        
        k = k+1;
    end
end

mask = ones(size(xpos));
id = find(isnan(xpos));
mask(id) = 0;

figure,

subplot(1,2,1)
plotretmap(xpos)
title('Horizontal retinotopy')
subplot(1,2,2)
plotretmap(ypos)
title('Vertical retinotopy')



%%
xposdum = xpos(2:end-1,2:end-1);
mi = prctile(xposdum(:),2); ma = prctile(xposdum(:),98);
id = find(xposdum<mi); xposdum(id) = mi;
id = find(xposdum>ma); xposdum(id) = ma;

yposdum = ypos(2:end-1,2:end-1);
mi = prctile(yposdum(:),2); ma = prctile(yposdum(:),98);
id = find(yposdum<mi); yposdum(id) = mi;
id = find(yposdum>ma); yposdum(id) = ma;

figure
subplot(1,3,1)
plotretmap(xposdum)
title('Horizontal retinotopy')
axis off

subplot(1,3,2)
plotretmap(yposdum)
title('Vertical retinotopy')
axis off

% subplot(1,3,3)
% contour(xposdum,'r'), axis square
% hold on
% contour(yposdum,'k'), axis square
% axis ij
% axis off


subplot(1,3,3)
plotorimap(omap), axis image

%Create the orientation legend%%%%%%%%%%%%%%%%%%
legdom = 0:30:180;
hsvdom = hsv;
id = round(linspace(1,64,length(legdom)));
hsvdom = hsvdom(id,:);
R = 1;
rid = linspace(1,length(xposdum(:,1)),length(legdom));
cid = length(xposdum(1,:)) + 5;
xpts_o = [0 0];
ypts_o = [1-R 1+R];

for i = 1:length(legdom)
   
    xpts = xpts_o*cos(legdom(i)*pi/180) + ypts_o*sin(legdom(i)*pi/180);
    ypts = xpts_o*sin(legdom(i)*pi/180) - ypts_o*cos(legdom(i)*pi/180);
    ypts = ypts + rid(i);
    xpts = xpts + cid;
    hold on
    line(xpts,ypts,'Color',hsvdom(i,:),'Clipping','off','LineWidth',3);
    
end
hold off

axis off
title('Orientation')
%work in progress...

% dori = 180/dim(1);
% oridom = 0:dori:180-dori;
% orimap = zeros(dim2(1),dim2(2));
% 
% 
% mag = abs(orimap);
% ang = angle(orimap)*180/pi;
% ang = (ang + (1-sign(ang))*180)/2;
% 
% %This is because the edges are wierd...
% mag = mag(3:end-2,3:end-2); ang = ang(3:end-2,3:end-2);
% 
% mag = mag-prctile(mag(:),0);
% mag = mag/prctile(mag(:),100);
% figure, imagesc(ang,'AlphaData',mag), colormap hsv
% 
% %Plot the map with anatomy
% figure
% dim = size(ang);
% set(gcf,'Color',[1 1 1]);
% 
% anatflag = 0;
% 
% if anatflag
%     
%     imanat = maskS.im{1};
%     imanat = imanat(3:end-2,3:end-2);
%    
%     mi = prctile(imanat(:),0);    
%     imanat = phi(imanat-mi);
%     ma = prctile(imanat(:),100);
%     id = find(imanat>ma);
%     imanat(id) = ma;
%     imanat = imanat/ma;
%     
%     mag = sqrt(imanat.*mag);
% 
%     imfunc = ang;
%     imfunc = imfunc/180;
%     imfunc = round(imfunc*63+1);
%     %imanat = round(imanat*63+1);
% 
%     hsvid = hsv;
%     imout = zeros(dim(1),dim(2),3);
%     for i = 1:dim(1)
%         for j = 1:dim(2)            
%             imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
%         end
%     end
%     
%     imanat(:,:,2) = imanat;
%     imanat(:,:,3) = imanat(:,:,1);
%     
%     imout = imout+(imanat).^.3;
% 
%     imout = imout/max(imout(:));
%     
%     %imout = imout(5:end,1:end-7,:);
%     
%     x = image(imout,'CDataMapping','direct','AlphaDataMapping','none');
% 
% else
%     imout = ang;
%     imout = imout/180;
%     imout = round(imout*63+1);
%     x = image(1:length(ang(1,:)),1:length(ang(:,1)),imout,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
%     
% end
% 
% axis image;
% 
% fh = gcf;
% 
% colormap hsv
% %colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})
% 
% %Create the orientation legend%%%%%%%%%%%%%%%%%%
% legdom = 0:30:180;
% hsvdom = hsv;
% id = round(linspace(1,64,length(legdom)));
% hsvdom = hsvdom(id,:);
% R = 4;
% rid = linspace(1,length(imout(:,1,1)),length(legdom));
% cid = length(imout(1,:,1)) + 8;
% xpts_o = [0 0];
% ypts_o = [1-R 1+R];
% 
% for i = 1:length(legdom)
%    
%     xpts = xpts_o*cos(legdom(i)*pi/180) + ypts_o*sin(legdom(i)*pi/180);
%     ypts = xpts_o*sin(legdom(i)*pi/180) - ypts_o*cos(legdom(i)*pi/180);
%     ypts = ypts + rid(i);
%     xpts = xpts + cid;
%     hold on
%     line(xpts,ypts,'Color',hsvdom(i,:),'Clipping','off','LineWidth',3);
%     
% end
% hold off
% 
% %%%%%%%%%%%%%%%


function plotretmap(im)

dom = linspace(min(im(:)),max(im(:)),5);
dom = round(dom*100)/100;

imfunc = im;
imfunc = imfunc-min(imfunc(:));
imfunc = imfunc/max(imfunc(:));
imfunc = round(imfunc*63+1);
dim = size(im);

jetid = jet;
imout = zeros(dim(1),dim(2),3);
for i = 1:dim(1)
    for j = 1:dim(2)
        if ~isnan(imfunc(i,j))
            imout(i,j,:) = jetid(imfunc(i,j),:);
        else
            imout(i,j,:) = [1 1 1];
        end
    end
end

image(imout)
axis image


for i = 1:length(dom)
    domcell{i} = dom(i);
end
iddom = linspace(1,64,length(dom));
colorbar('YTick',iddom,'YTickLabel',domcell)



function plotorimap(im)

dom = 0:45:180;

imfunc = im;
imfunc = imfunc/180;
imfunc = round(imfunc*63+1);
dim = size(im);

hsvid = hsv;
imout = zeros(dim(1),dim(2),3);
for i = 1:dim(1)
    for j = 1:dim(2)
        if ~isnan(imfunc(i,j))
            imout(i,j,:) = hsvid(imfunc(i,j),:);
        else
            imout(i,j,:) = [1 1 1];
        end
    end
end

image(imout)
axis image

