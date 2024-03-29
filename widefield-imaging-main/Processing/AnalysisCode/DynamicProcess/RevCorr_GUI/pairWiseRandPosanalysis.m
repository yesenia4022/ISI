function pairWiseRandPosanalysis

%Ian Nauhaus

global TC MK PW

global posprincax oriprincax

[xmicperpix ymicperpix] = getImResolution;

doriAll = []; dposAll = []; sizesumAll = []; doriNormAll = []; dposNormAll = []; DistAll = []; axAll = [];
for i = 1:MK.Ncell
    
    for j = i+1:MK.Ncell
        
        %Don't take abs() of dori/dpos... we want the sign to compute the
        %gradient direction
        dori = (oridiff(TC.OAng{1}(i)*pi/180,TC.OAng{1}(j)*pi/180)*180/pi); %degrees        
        dpos = sqrt((TC.xpos{1}(i)-TC.xpos{1}(j))^2 + (TC.ypos{1}(i)-TC.ypos{1}(j))^2);
        normer = sqrt(min([TC.xsize{1}(i) min(TC.ysize{1}(i))])*min([TC.xsize{1}(j) min(TC.ysize{1}(j))]));
        dpos = dpos/normer;

        doriAll = [doriAll dori];
        dposAll = [dposAll dpos];
        
        doriNorm = 2*dori/(TC.OSig{1}(i) + TC.OSig{1}(j));
        
        size1 = sqrt(TC.xsize{1}(i)*TC.ysize{1}(i));
        size2 = sqrt(TC.xsize{1}(j)*TC.ysize{1}(j));
        %size1 = min([TC.xsize{1}(i) TC.ysize{1}(i)]); 
        %size2 = min([TC.xsize{1}(j) TC.ysize{1}(j)]);
        
        sizesum = size1 + size2;
        dposNorm = 2*dpos/sizesum;
        
        doriNormAll = [doriNormAll doriNorm];
        dposNormAll = [dposNormAll dposNorm];
        
        sizesumAll = [sizesumAll sizesum];
        
        dy = (MK.CoM(i,1)-MK.CoM(j,1))*ymicperpix;
        dx = (MK.CoM(i,2)-MK.CoM(j,2))*xmicperpix;
        Dist = sqrt(dy^2 + dx^2); %Dist between cells in microns
        DistAll = [DistAll Dist];
        
        axAll = [axAll atan2(dy,dx)*180/pi];  %we want atan2 to compute gradient direction
        
    end
end
%doriAll = sizesumAll;

%%

Ddom = 0:65:250;
%Ddom = [0 40 92 160 252]; %logarithmic spacing
clear dori dpos doriNorm dposNorm sizesum
for i = 1:length(Ddom)-1

    DLimitMin = Ddom(i);  %Limit pairs to be this distance (microns) apart
    DLimitMax = Ddom(i+1);  %Limit pairs to be this distance (microns) apart
    id = find(DistAll<DLimitMax & DistAll>DLimitMin & ~isnan(doriAll.*dposAll) & ~isinf(doriAll.*dposAll));

    dori{i} = doriAll(id);
    dpos{i} = dposAll(id);

    doriNorm{i} = doriNormAll(id);
    dposNorm{i} = dposNormAll(id);
    
    sizesum{i} = sizesumAll(id);
    
    ax{i} = axAll(id);
    
    dist{i} = DistAll(id);

end

%% Get preferred map axis
%I already do this for all ROIs in posOriaxes.  I just put it here so that I
%could see it with the maps
axAll = [axAll axAll+180];
axAll = angle(exp(1i*axAll*pi/180))*180/pi;
doriAll = [doriAll -doriAll];
dposAll = [dposAll -dposAll];
DistAll = [DistAll DistAll];

id = find(DistAll<100);
DistAll = DistAll(id); axAll = axAll(id); doriAll = doriAll(id); dposAll = dposAll(id);

axedges = -180:20:180; dax = axedges(2)-axedges(1);
axW = 2*dax;
axdom = axedges(1:end-1)+dax/2;
doriAlln = doriAll./DistAll*1000; %deg/mm
dposAlln = dposAll./DistAll*1000; %oct/mm
%      doriAlln = doriAll;
%      dposAlln = dposAll;

varflag = 1;  %if we want it to be non-directional (axis), set to 1
f = 1;
if varflag
    doriAlln = abs(doriAlln);
    dposAlln = abs(dposAlln);
    f = 2;
end

for j = 1:length(axdom)  %loop through axis domain

    id = find(axAll > axdom(j)-axW/2 & axAll < axdom(j)+axW/2 );

    dori_mu(j) = trimmean(doriAlln(id),30);
    dori_SE(j) = std(doriAlln(id))/sqrt(length(id));

    dpos_mu(j) = trimmean(dposAlln(id),30);
    dpos_SE(j) = std(dposAlln(id))/sqrt(length(id));

end

% [mat xdom ydom] = smoothscatter(axAll,doriAlln,8,20,[-180 180],[prctile(doriAlln,2) prctile(doriAlln,98)]);
% [dum id] = max(mat);
% dori_mu = ydom(id);
% 
% id = find(dposAlln~=0);
% [mat xdom ydom] = smoothscatter(axAll(id),dposAlln(id),4,1,[-180 180],[prctile(dposAlln,2) prctile(dposAlln,98)]);
% [dum id] = max(mat);
% dpos_mu = ydom(id);
% axdom = xdom;


Resori = sum(dori_mu.*exp(1i*axdom*pi/180*f))/(.5*length(dori_mu));
oriprincax = angle(Resori)*180/pi/f;
Respos = sum(dpos_mu.*exp(1i*axdom*pi/180*f))/(.5*length(dpos_mu));
posprincax = angle(Respos)*180/pi/f;
M = abs(Respos.*Resori);
daxis = abs(oridiff(posprincax*pi/180,oriprincax*pi/180)*180/pi);

figure,
subplot(2,1,1)
scatter(axAll,doriAlln,'.')
hold on
%errorbar(axdom,dori_mu,dori_SE,'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
plot(axdom,dori_mu,'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
title(['phase = ' num2str(oriprincax) ' amp = ' num2str(abs(Resori))])
ylabel('deg/mm')

subplot(2,1,2)
scatter(axAll,dposAlln,'.')
hold on
%errorbar(axdom,dpos_mu,dpos_SE,'k','LineWidth',3), ylim([prctile(dposAlln,2) prctile(dposAlln,98)])
plot(axdom,dpos_mu,'k','LineWidth',3), ylim([prctile(dposAlln,2) prctile(dposAlln,98)])
title(num2str(posprincax))
title(['phase = ' num2str(posprincax) ' amp = ' num2str(abs(Respos)) '  diff = ' num2str(daxis)])
ylabel('oct/mm')
%%
[mat xdom ydom] = smoothscatter(abs(dori{1}),abs(dpos{1}),.8,.05);


figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori'), ylabel('dpos')
subplot(1,2,1),scatter(abs(dori{1}),abs(dpos{1}),'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori'), ylabel('dpos')
[r p] = corrcoef(abs(dori{1}),abs(dpos{1}));
title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])

%%
[mat xdom ydom] = smoothscatter(abs(doriNorm{1}),dposNorm{1},.015,.015);

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori (Norm dist)'), ylabel('dpos (Norm dist)')
subplot(1,2,1),scatter(abs(doriNorm{1}),dposNorm{1},'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori (Norm dist)'), ylabel('dpos (Norm dist)')
[r p] = corrcoef(abs(doriNorm{1}),abs(dposNorm{1}));
title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])

PW.dori = dori;
PW.dpos = dpos;
PW.sizesum = sizesum;
PW.Ddom = Ddom;
PW.doriNorm = doriNorm;
PW.dposNorm = dposNorm;
PW.ax = ax;
PW.dist = dist;

%%
Did = 1;
[mat xdom ydom] = smoothscatter(sizesum{Did},abs(dpos{Did}),.01,.01);

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('size sum'), ylabel('dpos')
subplot(1,2,1),scatter(sizesum{Did},abs(dpos{Did}),'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('size sum'), ylabel('dpos')
[r p] = corrcoef(sizesum{Did},abs(dpos{Did}));
title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])

%%
Did = 1;
[mat xdom ydom] = smoothscatter(sizesum{Did},abs(dori{Did}),.01,.01);

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('size sum'), ylabel('dori')
subplot(1,2,1),scatter(sizesum{Did},abs(dori{Did}),'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('size sum'), ylabel('dori')
[r p] = corrcoef(sizesum{Did},abs(dori{Did}));
title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])

function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

