function [maporidiff mapslopemag] = plotRCmaps_color

%Ian Nauhaus

global TC DM MK ACQinfo maskS idExamp

[xmicperpix ymicperpix] = getImResolution;

%Now plot color image of tuning
colorSensIM = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
phaseDiffIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
LMprojIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);

%Resample images to have equal resolution on both axes
xdom = (0:ACQinfo.pixelsPerLine-1)*xmicperpix;
ydom = (0:ACQinfo.linesPerFrame-1)*ymicperpix;
%%

xysize = [.33*xdom(end)/200 .34*ydom(end)/250];  %Size of the panels, normalized by the ROI size

c = 1;

for p = 1:MK.Ncell

    idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));
    
    colorSensIM(idcell) = TC.lumpref{1}(p)/2;
    
    phaseDiffIm(idcell) = TC.LMphaseDiff{1}(p);
    
    LMprojIm(idcell) = TC.LMproj{1}(p);

end


%%
%S vs. L+M
mag = abs(sign(colorSensIM));  id = find(isnan(mag)); mag(id) = -1;
IMtens = getImTens(colorSensIM,mag,[-.5 .5],'jet',1);
figure, subplot(1,3,1), image(xdom,ydom,phi(IMtens)), title('S-L-M/S+L+M'), axis image
plotEXcirc(idExamp,xdom,ydom)
F1F0dom = {-.5 ,-.25 ,0 ,.25, .5};
iddom = linspace(1,64,length(F1F0dom));
colorbar('YTick',iddom,'YTickLabel',F1F0dom)

%L vs. M phase difference
mag = sign(phaseDiffIm);  id = find(isnan(mag)); mag(id) = 0;
IMtens = getImTens(phaseDiffIm,mag,[0 180],'jet',1);  %this does not wrap, so use jet 
subplot(1,3,2), image(xdom,ydom,phi(IMtens)), title('L/M spatial phase difference'), axis image
plotEXcirc(idExamp,xdom,ydom)
phasedom = {0 ,45 ,90 ,135, 180};
iddom = linspace(1,64,length(phasedom));
colorbar('YTick',iddom,'YTickLabel',phasedom)

%L vs. M projection
mag = abs(sign(LMprojIm));  id = find(isnan(mag)); mag(id) = 0;
IMtens = getImTens(LMprojIm,mag,[-1 1],'jet',1);
subplot(1,3,3), image(xdom,ydom,phi(IMtens)), title('LM projection of phase vectors'), axis image
plotEXcirc(idExamp,xdom,ydom)
phasedom = {-1 ,0 ,1};
iddom = linspace(1,64,length(phasedom));
colorbar('YTick',iddom,'YTickLabel',phasedom)



function plotEXcirc(idExamp,xdom,ydom)

global MK

if ~isempty(idExamp)
    for q = 1:length(idExamp)
        hold on
        plot(xdom(round(MK.CoM(idExamp(q),2))),ydom(round(MK.CoM(idExamp(q),1))),'ok','MarkerSize',10)
    end
end

function IMtens = getImTens(pref,mag,mima,maptype,BackG)

id = find(mag>1);
mag(id) = 1;

id = find(pref<mima(1));
pref(id) = mima(1);
id = find(pref>mima(2));
pref(id) = mima(2);


prefid = (pref-mima(1))/(mima(2)-mima(1));
prefid = round(prefid*63+1);  %normalize to be colormap index

dim = size(pref);
if strcmp('blue',maptype)
    mapvals = gray;
    mapvals(:,1:2) = 0;
else
    mapvals = eval(maptype);
end

IMtens = zeros(dim(1),dim(2),3);
for i = 1:dim(1)
    for j = 1:dim(2)
        if mag(i,j) == 0 | isnan(mag(i,j)) | isnan(prefid(i,j));
            IMtens(i,j,:) = [BackG BackG BackG];
        else
            IMtens(i,j,:) = mag(i,j)*mapvals(prefid(i,j),:);
        end
    end
end
