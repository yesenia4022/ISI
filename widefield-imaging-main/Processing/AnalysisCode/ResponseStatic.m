function [ang] = ResponseStatic(Ncycles,f0dum)

global bw
 
lp = fspecial('gaussian',size(f0dum{5}),5);
lp = lp./sum(lp(:));
hp = zeros(size(lp));
rad = 50;
dum = fspecial('disk',rad);
hp(1:length(dum(:,1)),1:length(dum(1,:))) = dum;
hp = hp./sum(hp(:));

%HH = abs(fft2(lp))-abs(fft2(hp));
HH = abs(fft2(lp));
%HH = [];

[bwCell] = MakeCellMask;
bw = bw.*bwCell;
%bw = bwCell;

bwdum = double(bw);
id = find(bw(:) == 0);
bwdum(id) = NaN;

W = 5;

locR = [73 170];
locB = [119 232];

locR = [207 62]; %x y
locB =  [176 128]; %x y
xranR = (locR(1)-floor(W/2)):(locR(1)+floor(W/2));
yranR = (locR(2)-floor(W/2)):(locR(2)+floor(W/2));
xranB = (locB(1)-floor(W/2)):(locB(1)+floor(W/2));
yranB = (locB(2)-floor(W/2)):(locB(2)+floor(W/2));

k = 1;
for(i=0:length(f0dum)-1)
    pepsetcondition(i)
    if(~pepblank)       %This loop filters out the blanks
        v = pepgetvalues;
        oridomain(k) = v(1);
        
        id = find(isnan(f0dum{i+1}));
        f0dum{i+1}(id) = nanmedian(f0dum{i+1}(:));
        
        if ~isempty(HH)
            f0dum{i+1} = real(ifft2(HH.*fft2(f0dum{i+1})));    
        end
        
        tcB(k) = sum(sum(f0dum{i+1}(yranB,xranB)));
        tcR(k) = sum(sum(f0dum{i+1}(yranR,xranR)));

        f0{k} = Cellurize(f0dum{i+1},bwCell).*bw;
        

        %f0{k} = bwdum.*f0dum{i+1};
       
        k = k+1;
    end
end

%Put images in the right order
[oridomain id] = sort(oridomain);
tcB = tcB(id);
tcR = tcR(id);
for i = 1:length(id)
   f02{i} = f0{id(i)};
end
f0 = f02; 
clear f02

tcB = (tcB-min(tcB))/(max(tcB)-min(tcB));
tcR = (tcR-min(tcR))/(max(tcR)-min(tcR));

N = length(f0{1}(1,:));

%Normalize images to be 0 to 1
for i = 1:length(f0)
    ma(i) = max(f0{i}(:));
    mi(i) = min(f0{i}(:));
end
ma = max(ma);
mi = min(mi);
for i = 1:length(f0)
    f0{i} = (f0{i}-mi)/(ma-mi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%Normalize the tuning curve of each pixel

miIm = ones(N^2,1);
maIm = zeros(N^2,1);
for i = 1:length(f0)
    miIm = min([f0{i}(:) miIm(:)],[],2);
    maIm = max([f0{i}(:) maIm(:)],[],2);
end
maIm = reshape(maIm,N,N);
miIm = reshape(miIm,N,N);

for i = 1:length(f0)
    f0{i} = (f0{i}-miIm)./(maIm-miIm);
    %f0{i} = (f0{i}-miIm);
end

% sumIm = zeros(size(f0{1}));
% for i = 1:length(f0)
%     sumIm = sumIm + f0{i};
% end
% 
% for i = 1:length(f0)
%     f0{i} = f0{i}./sumIm;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Normalize images to be 0 to 1
for i = 1:length(f0)
    ma(i) = max(f0{i}(:));
    mi(i) = min(f0{i}(:));
end
ma = max(ma);
mi = min(mi);
for i = 1:length(f0)
    f0{i} = (f0{i}-mi)/(ma-mi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[x y] = meshgrid(0:N-1,0:N-1);

mid = ceil(N/2);
r = sqrt((x-mid).^2 + (y-mid).^2);
id = find(r<=mid*.75);
Mask = zeros(N,N);
Mask(id) = 1;


figure
plot(oridomain,tcB,'b','linewidth',2), hold on, plot(oridomain,tcB,'ob','linewidth',2)
hold on
plot(oridomain,tcR+1.3,'r','linewidth',2), hold on, plot(oridomain,tcR+1.3,'or','linewidth',2)
xlim([oridomain(1) oridomain(end)]), ylim([-.2 2.5]),
hold on
plot([oridomain(1) oridomain(end)],[1.15 1.15],'k')
title('Single Cell Tuning Curves','FontWeight','Bold'), xlabel('Orientation')
hold off

%%%%%%%%%%%%%%%%%%%%%%
%Make orientation map
%%%%%%%%%%%%%%%%%%%%%%
clear im

orimap = zeros(size(f0{1}));
for z = 1:length(f0)
    orimap = orimap + f0{z}*exp(1i*2*oridomain(z)*pi/180);    
end
ang = angle(orimap);
ang = (ang+pi*(1-sign(ang)))*180/pi/2;

id = find(isnan(ang));
ang(id) = 0;

cidx = hsv;
idx = 1+round(63*ang/180);
for i = 1:length(ang(:,1))
    for j = 1:length(ang(1,:))
        im(i,j,:) = squeeze(cidx(idx(i,j),:));
        if ~bw(i,j)
            im(i,j,:) = 1;
        end
    end
end

rrange = 50:150;
crange = 150:225;

figure,image(im(rrange,crange,:))
hold on
plot(locB(1)-crange(1)+1,locB(2)-rrange(1)+1,'bo','markersize',25,'linewidth',3)
hold on
plot(locR(1)-crange(1)+1,locR(2)-rrange(1)+1,'ro','markersize',25,'linewidth',3)
hold off

set(gca,'Xtick',[],'Ytick',[])
set(gcf,'Color',[1 1 1]);
title('Oriention Preferences','FontWeight','Bold')
MakeOriBars(ang(rrange,crange))
