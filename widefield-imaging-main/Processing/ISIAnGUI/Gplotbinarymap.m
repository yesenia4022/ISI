function Gplotbinarymap(binMap,anatflag)

global fh G_handles Analyzer bw
%mag = log(mag)

%mag = ones(size(mag));

dim = size(binMap);
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
        
        if isempty(Analyzer.ACQ.ROIcrop)
            yran = 1:size(im,1);
            xran = 1:size(im,2);
        else
            yran = Analyzer.ACQ.ROIcrop(2):(Analyzer.ACQ.ROIcrop(4)+Analyzer.ACQ.ROIcrop(2)-1);
            xran = Analyzer.ACQ.ROIcrop(1):(Analyzer.ACQ.ROIcrop(3)+Analyzer.ACQ.ROIcrop(1)-1);
        end
        
        im = im(yran,xran);
        imanat = double(im);
        
    end
    
    mi = prctile(imanat(find(bw)),10);
    imanat = phi(imanat-mi);
    ma = prctile(imanat(find(bw)),90);
    imanat = imanat/ma;
    %%
%     figure,imagesc(imfunc), colormap gray, 
%     bwV1 = roipoly;
%     
%     figure,imagesc(imanat,[0 2]), colormap gray, hold on, contour(bwV1,[.5 .5],'r')
%     axis square
%     axis off
%     figure,imagesc(imfunc), colormap gray, hold on, contour(bwV1,[.5 .5],'r')
%     axis square
%     axis off
    
    %%%
    mag = (imanat.*bw).^.5;
    %%%
    
    imfunc = binMap;
    
    mi = prctile(imfunc(find(bw)),10); 
    ma = prctile(imfunc(find(bw)),90);
    imfunc = (imfunc-mi) / (ma-mi);
    id = find(imfunc(:)<0);
    imfunc(id) = 0;
    id = find(imfunc(:)>1);
    imfunc(id) = 1;

    imfunc = round(imfunc*63+1);
    %imanat = round(imanat*63+1);
    
    
    
    hsvid = jet;
    imout = zeros(dim(1),dim(2),3);
    for i = 1:dim(1)
        for j = 1:dim(2)            
            imout(i,j,:) = mag(i,j)*hsvid(imfunc(i,j),:);
        end
    end
    
    imanatRGB = imanat;
    imanatRGB(:,:,2) = imanatRGB;
    imanatRGB(:,:,3) = imanatRGB(:,:,1);
    
    %imout = 3*imout.^3+.3*(imanat).^.3;
    
    imout = imout + imanatRGB;

    imout = imout/max(imout(:));
    
    %imout = imout(1:end-4,8:end,:);
    figure
    x = image(imout,'CDataMapping','direct','AlphaDataMapping','none');

else
    figure
    imout = binMap;
%     imout = imout/180;
%     imout = round(imout*63+1);
%bw(find(isnan(imout))) = 0;

    mi = prctile(imout(find(bw)),10); 
    ma = prctile(imout(find(bw)),90);

    imagesc(imout.*bw,'CDataMapping','direct','AlphaData',bw,'AlphaDataMapping','none',[mi ma]);
    colormap jet

    
end

axis image;

fh = gcf;
colorbar

colormap jet


%%%%%%%%%%%%%%%


datacursormode on;
dcm_obj = datacursormode(fh);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','UpdateFcn',@myupdatefcn);


function txt = myupdatefcn(empt,event_obj)

%Matlab doesn't like it when I try to input other things into myupdatefcn,
%this is why I have these globals
 
global Analyzer ACQinfo Tens Tens_var Flim TCWin Fsymbol G_handles oppCollapse

figure(99)

W = TCWin;

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

[tc tcourseHi tcourseLo tcourseHi_var tcourseLo_var tcourseBlank axisdom blank legStr] = getpixeldata2(pos,W);  %This does the work
%%%


subplot(2,1,1)
if ~isempty(blank)
    plot([axisdom(1) axisdom(end)],[blank blank],'k'), hold on
else 
    %Even if no blank was shown we put a line at zero. 
    plot([axisdom(1) axisdom(end)],[0 0],'k'), hold on  
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
% if ~isempty(blank)
%     tc = tc-blank;
% end

d = size(tc);
if d(1) == 1 || d(2) == 1
    tcdum = tc;
else
    [y x] = find(tc == max(tc(:)));
    tcdum = tc(:,x);
end


% 
% TCSel = abs(tc(:)'*exp(1i*2*axisdom(:)*pi/180));
% TCSel = TCSel./sum((tc))';
% TCSel =  round(TCSel*100)/100;
% title(['OSI = ' num2str(TCSel')])

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

%This is for if the shutter opens/closes at the beginning/end
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
       ['Ori: ' sprintf('%2.1f %%',data(round(pos(2)),round(pos(1)))/64*180) ' deg']};
       