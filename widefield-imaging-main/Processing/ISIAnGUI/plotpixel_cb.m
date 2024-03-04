function plotpixel_cb

fh = gcf;
%colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','45','90','135','180'})

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global ACQinfo Tens f0m Analyzer TCWin

W = TCWin;

tdom = 0:length(Tens{1}(1,1,:))-1;
tdom = tdom*ACQinfo.msPerLine/1000*ACQinfo.linesPerFrame;
predelay = getparam('predelay');
trialtime = getparam('stim_time');
tdom = tdom-predelay;

screenDist = Analyzer.M.screenDist;
screenResX = Analyzer.M.xpixels/Analyzer.M.screenXcm;  %pix/cm
screenResY = Analyzer.M.ypixels/Analyzer.M.screenYcm;

[xpos ypos xsize ysize] = getPosSize;
x_size = (max(xpos)-min(xpos)+min(xsize))/screenResX;  %cm stimulus width
x_size = 2*atan2(x_size/2,screenDist)*180/pi;  %convert to deg
y_size = (max(ypos)-min(ypos)+min(ysize))/screenResY;  %cm stimulus width
y_size = 2*atan2(y_size/2,screenDist)*180/pi;  %convert to deg

rows = ACQinfo.linesPerFrame;
cols = ACQinfo.pixelsPerLine;
%  
pos = get(event_obj,'Position'); %pos(1) is column dimension (top left is origin)

xran = (pos(1)-floor(W/2)):(pos(1)+floor(W/2));
yran = (pos(2)-floor(W/2)):(pos(2)+floor(W/2));

tau = pos(2)*ACQinfo.msPerLine/1000;
tdom = tdom + tau;

figure(99)

nc = getnoconditions;
bflag = stimblank(nc);
if bflag
    nloop = nc-1;
else
    nloop = nc;
end
for j = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp('b',Analyzer.loops.conds{1}.symbol{j})
        bid = j;
    end
end

retorixy = [0 90];
for z = 1:2 %loop through x and y dimension
    clear phasedom
    tc = NaN*ones(1,nloop);
    for i = 1:nloop  %loop through each condition
        b = Analyzer.loops.conds{i}.val{bid};
        ori = round(90*b);

        if retorixy(z) == ori

            dum = f0m{i}(yran,xran);
            tc(i) = mean(dum(:));

        end
    end
    
    if retorixy(z) == 0
        phasedom = xpos;
    else
        phasedom = ypos;
    end
    
    phasedom(isnan(phasedom)) = [];

    [ma idma] = max(tc);
    [mi idmi] = min(tc);
    tc(find(isnan(tc))) = []; %Get rid of blank index

    [phasedom id] = sort(phasedom);

    subplot(2,2,z)
    plot(phasedom,tc(id))
    if z == 1         
        xlabel('x position')
    else        
        xlabel('y position')
    end
    
    nopix = length(yran)*length(xran);

    subplot(2,2,z+2)
    dum = squeeze(sum(sum(Tens{idma}(yran,xran,:),1),2))/nopix;
    plot(tdom(1:end-1),dum(1:end-1)), hold on, plot(tdom(1:end-1),dum(1:end-1),'o')
    hold on
    dum = squeeze(sum(sum(Tens{idmi}(yran,xran,:),1),2))/nopix;
    plot(tdom(1:end-1),dum(1:end-1),'r'), hold on, plot(tdom(1:end-1),dum(1:end-1),'or')
    ylimits = get(gca,'Ylim');
    hold on, plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')

    hold off
    xlabel('sec')


    tar = get(get(event_obj,'Target'));
    data = tar.CData;

    txt = {['X: ',num2str(pos(1))],...
        ['Y: ',num2str(pos(2))],...
        ['Pos: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180)]};
end