function Gplotdirmap(mag,ang,anatflag)

global fh
%mag = log(mag)
mag = mag-min(mag(:));
mag = mag/max(mag(:));

dim = size(ang);
set(gcf,'Color',[1 1 1]);

if anatflag
    
        
    CH = GetTrialData([-inf inf],1);
%     if get(G_handles.fastMotionFlag,'Value')
%         [Px_fast Py_fast] = getTrialMotion3(CH{1});
%         CH{1} = makeGeoTrx(CH{1},Px_fast,Py_fast);
%     end
    imanat = mean(CH(:,:,2:end-1),3);
    

    %% get anatomy pic
    [filename, pathname] = uigetfile([ '/*.mat'], 'Select a grab manually');
    
    if filename~=0
        
        S = load(strcat(pathname,filename));  %Returns the contents in the .mat under the structure S
        
        if isfield(S,'grab')
            im = S.grab.img;
        else
            im = S.im;
        end
        
        imanat = double(im);
        
    end
    
    mi = prctile(imanat(:),0);
    imanat = phi(imanat-mi);
    ma = prctile(imanat(:),100);
    imanat = imanat/ma;
    %%
    
    %%%
    mag = (imanat.*mag).^.5;
    %%%
    

    imfunc = ang;
    imfunc = imfunc/360;
    imfunc = round(imfunc*63+1);
    %imanat = round(imanat*63+1);

    hsvid = hsv;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
        end
    end
    
%     imanat(:,:,2) = imanat;
%     imanat(:,:,3) = imanat(:,:,1);
%     
%     imout = imout+sqrt(imanat);
% 
%     imout = imout/max(imout(:));

    figure
    
    x = image(imout,'CDataMapping','direct','AlphaDataMapping','none');

else
    figure
    
    imfunc = 64*ang/360;
    x = image(1:length(ang(1,:)),1:length(ang(:,1)),imfunc,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    colormap hsv
    
end

axis image;

fh = gcf;

colormap hsv
colorbar('YTick',[1 16:16:64],'YTickLabel',{'0','90','180','270','360'})

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global Analyzer ACQinfo Tens Tens_var Flim TCWin Fsymbol G_handles oppCollapse

W = TCWin;

figure(99)

varflag = get(G_handles.EbarFlag,'Value');
    
tdom = Analyzer.syncInfo{1}.acqSyncs;
tdom = tdom-tdom(1);

predelay = getParamVal('predelay',0);
trialtime = getParamVal('stim_time',0);
tdom = tdom-predelay;
nr = getnorepeats(1);

SEn = sqrt(length(Flim(1):Flim(2))*nr);  %standard error normalizer for tuning curve
%  
pos = round(get(event_obj,'Position')); %pos(1) is column dimension

%%%
[tc tcourseHi tcourseLo tcourseHi_var tcourseLo_var axisdom blank legStr] = getpixeldata(pos,W);  %This does the work
%%%


subplot(2,1,1)
if ~isempty(blank)
    plot([axisdom(1) axisdom(end)],[blank blank],'k'), hold on
else 
    %Even if no blank was shown we put a line at zero. 
    %plot([axisdom(1) axisdom(end)],[0 0],'k'), hold on  
end


plot(axisdom,tc,'-o'), hold off
legend(legStr)
% if ~varflag
%     plot(axisdom,tc,'-o'), hold off
%     legend(legStr)
% else
%     errorbar(axisdom,tc(id),sqrt(tc_var(id))/SEn,'b'), hold off
% end
xlabel(Fsymbol)


%Get 'orientation selectivity index' and put into the title
if ~isempty(blank)
    tc = tc-blank;
end

d = size(tc);
if d(1) == 1 || d(2) == 1
    tcdum = tc;
else
    [y x] = find(tc == max(tc(:)));
    tcdum = tc(:,x);
end

TCSel = abs(sum(tcdum'.*exp(1i*axisdom*pi/180)));
TCSel = TCSel/sum(tcdum);
TCSel =  round(TCSel*100)/100;
title(['DSI = ' num2str(TCSel)])


Fi = 1;

subplot(2,1,2)
if varflag
    %dum_var = squeeze(sum(sum(Tens_var{idma}(yran,xran,:),1),2))/nopix/nr;
    
    dumSE = sqrt(tcourseHi_var(1:end-Fi)/nr);
    errorbar(tdom(1:end-Fi),tcourseHi(1:end-Fi),dumSE), hold on 
else
    plot(tdom(1:end-Fi),tcourseHi(1:end-Fi),'.-'), hold on
end

if varflag
    %dum_var = squeeze(sum(sum(Tens_var{idmi}(yran,xran,:),1),2))/nopix/nr;
    
    dumSE = sqrt(tcourseLo_var(1:end-Fi)/nr);
    errorbar(tdom(1:end-Fi),tcourseLo(1:end-Fi),dumSE,'r')
else
    plot(tdom(1:end-Fi),tcourseLo(1:end-Fi),'.-r')
end

if isfield(ACQinfo,'stimPredelay')
    ylimits = get(gca,'Ylim');
    plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
end
hold off
xlabel('sec')



tar = get(get(event_obj,'Target'));
data = tar.CData;

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Dir: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*360) ' deg']};
       
       