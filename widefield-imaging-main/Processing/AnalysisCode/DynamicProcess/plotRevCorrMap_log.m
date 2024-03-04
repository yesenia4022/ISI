function plotRevCorrMap_log

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

t1 = 150; t2 = 500;
[dum id1] = min(abs(taudom-t1));
[dum id2] = min(abs(taudom-t2));

Tens = zeros(dim2(1),dim2(2),dim(2));
for ori = 1:dim(1)    
    for sf = 1:dim(2)
        for phase = 1:dim(3)
            Tens(:,:,sf) = Tens(:,:,sf) + mean(kernelsIm{ori,sf,phase}(:,:,id1:id2),3);
        end
    end    
end

h = fspecial('gaussian', size(Tens(:,:,1)), .5);
h = abs(fft2(h));
for i = 1:dim(2)   
    dum = Tens(:,:,i);
    
    id = find(isnan(dum(:)));
    dum(id) = 0;

    dum = ifft2(fft2(dum).*h);
    
    Tens(:,:,i) = dum;

end

[ma id] = max(Tens,[],3);
[mi dum] = min(Tens,[],3);
mag = (ma-mi);
mag = mag-min(mag(:));
mag = mag/max(mag(:));
figure,imagesc(id,'alphadata',mag)
axis image

mi = min(Tens,[],3);
for sf = 1:dim(2)   
    Tens(:,:,sf) = Tens(:,:,sf)-mi;
end
su = sum(Tens,3);
for sf = 1:dim(2)   
    Tens(:,:,sf) = Tens(:,:,sf)./su;
end

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