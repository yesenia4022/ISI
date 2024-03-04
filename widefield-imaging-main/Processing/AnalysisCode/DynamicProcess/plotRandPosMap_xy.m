function plotRandPosMap_xy

%Ian Nauhaus

global kernelsIm G_RChandles maskS ACQinfo

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

%Make a spatial tuning curve for each orientation
for ori = 1:dim(1)
    for bw = 1:2
        Tensbw{ori,bw} = zeros(dim2(1),dim2(2),dim(2));
        for x = 1:dim(2)
            Tensbw{ori,bw}(:,:,x) = Tensbw{ori,bw}(:,:,x) + mean(kernelsIm{ori,x,1,bw}(:,:,id1:id2),3);
        end
    end
    Tens{ori} = (Tensbw{ori,1} + Tensbw{ori,2})/2;
end

DS = 0;  %downsample
if DS
    for ori = 1:dim(1)
        Tens{ori} = (Tens{ori}(:,:,1:2:end) + Tens{ori}(:,:,2:2:end))/2;
    end
    xdom = xdom(2:2:end); xdom = xdom-mean(xdom);
    dim(2) = length(xdom);
    
    for i = 1:dim(1)/2
        Tens2{i} = Tens{i*2-1} + Tens{i*2};
    end
    Tens = Tens2;  clear Tens2
    oridom = (oridom(1:2:end) + oridom(2:2:end))/2;
    dim(1) = length(oridom);
    
end

h = fspecial('gaussian', size(Tens{1}(:,:,1)), 1);
h = abs(fft2(h));

k = 1;
clear v
%figure
for ori = 1:dim(1)
    ori
    for i = 1:dim(2)
        dum = Tens{ori}(:,:,i);

        dum = ifft2(fft2(dum).*h);
        
        Tens{ori}(:,:,i) = dum;
    end

    [ma id] = max(Tens{ori},[],3);
    [mi dum] = min(Tens{ori},[],3);
    mag = (ma-mi);
    mag = mag-min(mag(:));
    mag = mag/max(mag(:));

%     mi = median(Tens{ori},3);
%     for x = 1:dim(2)
%         Tens{ori}(:,:,x) = phi(Tens{ori}(:,:,x)-mi);
%     end
    
    oritc(:,:,ori) = std(Tens{ori},[],3); %effectively the ori tuning curve at each pixel
    
%     su = sum(Tens{ori},3);
%     for x = 1:dim(2)
%         Tens{ori}(:,:,x) = Tens{ori}(:,:,x)./su;
%     end
    
    for x = 1:dim(2)
        %subplot(dim(1),dim(2),k)
        %imagesc(Tens{ori}(:,:,x),[prctile(Tens{ori}(:),1) prctile(Tens{ori}(:),99)])
        %plot(taudom,squeeze(mean(mean(kernelsIm{k},1),2))), ylim([0 20])
        
        blockSize = 8;
        Nblocks = ACQinfo.linesPerFrame/blockSize;
        for q = 1:Nblocks  %get position/orientation kernel for each block
            for p = 1:Nblocks
                dum = Tens{ori}(1+(q-1)*blockSize:q*blockSize,1+(p-1)*blockSize:p*blockSize,x);
                v{q,p}(x,ori) = mean(dum(:)); 
            end
        end
        
        
        drawnow
        k = k+1;
    end
    
    %We get an estimate of the x and y position for each orientation...
    Pref{ori} = 0;

    
end

su = sum(oritc,3);
for ori = 1:dim(1)
    oritc(:,:,ori) = oritc(:,:,ori)./su;
end

%%

degperpix = getparam('x_size')/(length(v(:,1))-1);  %can't use xdom because of the padding

clear RF xpos ypos
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
        RF{q,p}(:,1:end/2) = median(RF{q,p}(:));
        
%         xdom = (0:length(RF(1,:))-1)*degperpix; %degrees
%         ydom = (0:length(RF(:,1))-1)*degperpix;   %don't use 'y_size'

        thresh2 = prctile(RF{q,p}(:),20) + (max(RF{q,p}(:))-prctile(RF{q,p}(:),20))*.3;

        RF{q,p} = phi(RF{q,p}-thresh2);

        otc = max(vdum);

        [dum omap(q,p)] = orifind(otc',oridom);

        [param ffit varaccount] = Gaussfit2Drot(RF{q,p},omap(q,p));

        xpos(q,p) = param(2)*degperpix;
        ypos(q,p) = param(1)*degperpix;

%         subplot(Nblocks,Nblocks,k)
%         imagesc(RF{q,p}), drawnow
%         axis square
%         axis off
        
        k = k+1;
    end
end

figure,

subplot(1,2,1)
imagesc(xpos), axis square
colorbar
title('Horizontal retinotopy')
subplot(1,2,2)
imagesc(ypos), axis square
colorbar
title('Vertical retinotopy')



%u/v are brain coordinates; x/y are retinotopy
dxdu = xpos(:,2:end) - xpos(:,1:end-1);
dxdv = xpos(2:end,:) - xpos(1:end-1,:);
dydu = ypos(:,2:end) - ypos(:,1:end-1);
dydv = ypos(2:end,:) - ypos(1:end-1,:);

dposdu = sqrt(dxdu.^2 + dydu.^2);
dposdv = sqrt(dxdv.^2 + dydv.^2);

doriu = abs(oridiff(omap(:,2:end)*pi/180,omap(:,1:end-1)*pi/180));
doriv = abs(oridiff(omap(2:end,:)*pi/180,omap(1:end-1,:)*pi/180));

figure,
subplot(1,2,1),scatter(dposdu(:),doriu(:),'.')
[r p] = corrcoef(dposdu(:),doriu(:));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])
subplot(1,2,2),scatter(dposdv(:),doriv(:),'.')
[r p] = corrcoef(dposdv(:),doriv(:));
title(['r = ' num2str(r(1,2))  ';  p = ' num2str(p(1,2))])

%Compute x and y position based on a weighted average the estimates at each
%orientation.  The weighing function is effectively the ori tuning curve.
%%
% xpos = 0;
% ypos = 0;
% 
% id = find(abs(cos(oridom*pi/180))>=.7);
% for i = 1:length(id)   
%     xpos = xpos + cos(oridom(id(i))*pi/180)*Pref{id(i)}/length(id);
% end
% 
% id = find(abs(sin(oridom*pi/180))>=.7);
% for i = 1:length(id)   
%     ypos = ypos + sin(oridom(id(i))*pi/180)*Pref{id(i)}/length(id);
% end


% for ori = 1:dim(1)   
% %     xpos = xpos + cos(oridom(ori)*pi/180)*Pref{ori}.*oritc(:,:,ori);
% %     ypos = ypos + sin(oridom(ori)*pi/180)*Pref{ori}.*oritc(:,:,ori); 
%     xpos = xpos + cos(oridom(ori)*pi/180)*Pref{ori}/dim(1);
%     ypos = ypos + sin(oridom(ori)*pi/180)*Pref{ori}/dim(1);  
% end

xpos = xpos(2:end-1,2:end-1);
ypos = ypos(2:end-1,2:end-1);

figure,
subplot(1,3,1)
imagesc(xpos,[prctile(xpos(:),2) prctile(xpos(:),98)]), colorbar, axis square
title('xpos')
subplot(1,3,2)
imagesc(ypos,[prctile(ypos(:),2) prctile(ypos(:),98)]), colorbar, axis square
title('ypos')
subplot(1,3,3)
contour(xpos,'r'), axis square
hold on
contour(ypos,'k'), axis square
axis ij

figure,
imagesc(omap,[0 180]), colormap hsv, axis square, colorbar

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