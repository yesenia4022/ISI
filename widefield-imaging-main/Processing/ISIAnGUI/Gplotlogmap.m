function Gplotlogmap(mag,pref,anatflag,varargin)

global fh symbolInfo


if ~isempty(varargin)
    transparentMask = varargin{1};
    transparentID = find(transparentMask == 1);
end

logdom = getdomain(symbolInfo.str{1});

mag = mag-prctile(mag(:),1);
mag(find(mag<0)) = 0;
mag = mag/prctile(mag(:),99.5);
mag(find(mag>1)) = 1;

%mag = ones(size(mag));


pref = log10(pref);
dim = size(mag);
set(gcf,'Color',[1 1 1]);

id = find(isnan(pref));
mag(id) = 0;
pref(id) = min(pref(:));

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

    imfunc = pref;
    if logdom(1) == 0
        logdom(1) = logdom(2)/10;
    end
    imfunc = imfunc-log10(logdom(1));
    imfunc = imfunc/(log10(logdom(end))-log10(logdom(1)));
    imfunc = round(imfunc*63+1);  %domain max at 64
    
    jetid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)
            imout(i,j,:) = mag(i,j)*jetid(imfunc(i,j),:);
        end
    end
    
    imanat(:,:,2) = imanat;
    imanat(:,:,3) = imanat(:,:,1);
    
    imout = imout+sqrt(imanat);

    imout = imout/max(imout(:));
    figure
    x = image(imout,'CDataMapping','direct','AlphaDataMapping','none');
    
    
    for i = 1:length(logdom)
        domcell{i} = logdom(i);
    end
    iddom = linspace(0,1,length(logdom));
    colorbar('YTick',iddom,'YTickLabel',domcell)

else      

    imfunc = pref;
    if logdom(1) == 0
        logdom(1) = logdom(2)/10;
    end
    imfunc = imfunc-log10(logdom(1));
    imfunc = imfunc/(log10(logdom(end))-log10(logdom(1)));
    imfunc = round(imfunc*63+1);  %domain max at 64
    
    figure
    imagesc(imfunc,'CDataMapping','direct','AlphaData',mag,'AlphaDataMapping','none');
    colormap jet
    
    for i = 1:length(logdom)
        domcell{i} = logdom(i);
    end
    iddom = linspace(1,64,length(logdom));
    colorbar('YTick',iddom,'YTickLabel',domcell)
end


colormap jet
axis image;

fh = gcf;

datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global Analyzer ACQinfo Tens Tens_var Flim TCWin Fsymbol oppCollapse G_handles

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

%pos = [696 417]

%%%
[tc tcourseHi tcourseLo tcourseHi_var tcourseLo_var tcourseBlank logdom blank legStr] = getpixeldata2(pos,W);  %This does the work
%%%


subplot(2,1,1)
if ~isempty(blank)
    plot([logdom(1) logdom(end)],[blank blank],'k'), hold on
else 
    %Even if no blank was shown we put a line at zero. 
    plot([logdom(1) logdom(end)],[0 0],'k'), hold on  
end

plot(logdom,tc,'-o'), hold off
legend(legStr)

% if ~varflag
%     plot(logdom,tc,'-o'), hold off
%     legend(legStr)
% else
%     errorbar(logdom,tc(id),sqrt(tc_var(id))/SEn,'b'), hold off
% end
set(gca,'Xtick',round(sort(logdom)*1000)/1000,'Xscale','log')
id = find(Fsymbol ~= '_');
xlabel(Fsymbol(id))


%Get 'orientation selectivity index' and put into the title
if ~isempty(blank)
    tc = tc-blank;
end


if oppCollapse == 3
    [y x] = find(tc == max(tc(:)));
    tcdum = tc(:,x);
else
    tcdum = tc;
end


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
    hold on
    plot(tdom(1:end-Fi),tcourseBlank(1:end-Fi),'.-k')
end


dum = [tcourseLo(:); tcourseHi(:)];
ylim([prctile(dum,1)  prctile(dum,100)]);

ylimits = get(gca,'Ylim');
plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')

hold off
xlabel('sec')



tar = get(get(event_obj,'Target'));
data = tar.CData;

txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Value: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180)]};
       

       