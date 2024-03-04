function pairWiseRandPosanalysis2

%Ian Nauhaus

%2 doesn't parse cortical distance

global TC MK PW

global posprincax oriprincax

[xmicperpix ymicperpix] = getImResolution;

PW = struct;

Npair = 0;
for i = 1:MK.Ncell
    for j = (i+1):MK.Ncell
        Npair = Npair+1;
    end
end

doriAll = []; dposAll = []; sizesumAll = []; doriNormAll = []; dposNormAll = []; DistAll = []; axAll = [];

for c = 1:length(TC.xsize)  %length(DM.colordom) might be different
    PW.dori{c} = zeros(1,Npair); PW.dpos{c} = zeros(1,Npair); 
    PW.sizesum{c} = zeros(1,Npair); 
    PW.doriNorm{c} = zeros(1,Npair); PW.dposNorm{c} = zeros(1,Npair);    
    PW.Dist{c} = zeros(1,Npair); PW.ax{c} = zeros(1,Npair);

    k = 1;
    for i = 1:MK.Ncell

        for j = i+1:MK.Ncell

            %Don't take abs() of dori/dpos... we want the sign to compute the
            %gradient direction
            
            if (isfield(TC,'sfreq'))
                PW.dsf{c}(k) = abs(log2(TC.sfreq{1}(i)/TC.sfreq{1}(j)));
            end
            
            PW.dori{c}(k) = (oridiff(TC.OAng{1}(i)*pi/180,TC.OAng{1}(j)*pi/180)*180/pi); %degrees
            PW.dpos{c}(k) = sqrt((TC.xpos{1}(i)-TC.xpos{1}(j))^2 + (TC.ypos{1}(i)-TC.ypos{1}(j))^2);

            PW.doriNorm{c}(k) = 2*PW.dori{c}(k)/(TC.OSig{1}(i) + TC.OSig{1}(j));

%             size1 = sqrt(TC.xsize{c}(i)*TC.ysize{c}(i));
%             size2 = sqrt(TC.xsize{c}(j)*TC.ysize{c}(j));
            size1 = min([TC.xsize{c}(i) TC.ysize{c}(i)]);
            size2 = min([TC.xsize{c}(j) TC.ysize{c}(j)]);

            PW.sizesum{c}(k) = size1 + size2;
            PW.dposNorm{c}(k) = 2*PW.dpos{c}(k)/PW.sizesum{c}(k);
            
            PW.dsize{c}(k) = abs(log2(size1/size2));
            
            PW.dBWdiff{c}(k) = TC.BWdiff{c}(i) - TC.BWdiff{c}(j);

            dy = (MK.CoM(i,1)-MK.CoM(j,1))*ymicperpix;
            dx = (MK.CoM(i,2)-MK.CoM(j,2))*xmicperpix;
            PW.Dist{c}(k) = sqrt(dy^2 + dx^2); %Dist between cells in microns
            
            v1 = TC.RF{1}{i}(:); v2 = TC.RF{1}{j}(:);            
            dum = corrcoef(v1,v2); dum = -dum(1,2); PW.RFEuc{c}(k) = (dum+1)/2;

            PW.ax{c}(k) = atan2(dy,dx)*180/pi;
            k = k+1;
        end
    end

    %doriAll = sizesumAll;
    PW.dori{c} = abs(PW.dori{c});
    PW.doriNorm{c} = abs(PW.doriNorm{c});

%     PW.doriNorm{c} = (PW.doriNorm{c}-nanmean(PW.doriNorm{c}))/nanstd(PW.doriNorm{c});
%     PW.dori{c} = (PW.dori{c}-nanmean(PW.dori{c}))/nanstd(PW.dori{c});
%     PW.sizesum{c} = (PW.sizesum{c}-nanmean(PW.sizesum{c}))/nanstd(PW.sizesum{c});
%     PW.dpos{c} = (PW.dpos{c}-nanmean(PW.dpos{c}))/nanstd(PW.dpos{c});
%     PW.dposNorm{c} = (PW.dposNorm{c}-nanmean(PW.dposNorm{c}))/nanstd(PW.dposNorm{c});

end

%%
% [mat xdom ydom] = smoothscatter(abs(PW.dori{1}),abs(PW.dpos{1}),.8,.05);
% 
% 
% figure,
% subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
% xlabel('dori'), ylabel('dpos')
% subplot(1,2,1),scatter(abs(PW.dori{1}),abs(PW.dpos{1}),'.k')
% xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
% xlabel('dori'), ylabel('dpos')
% [r p] = corrcoef(abs(PW.dori{1}),abs(PW.dpos{1}));
% title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])
% 
% %%
% [mat xdom ydom] = smoothscatter(abs(PW.doriNorm{1}),PW.dposNorm{1},.015,.015);
% 
% figure,
% subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
% xlabel('dori (Norm dist)'), ylabel('dpos (Norm dist)')
% subplot(1,2,1),scatter(abs(PW.doriNorm{1}),PW.dposNorm{1},'.k')
% xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
% xlabel('dori (Norm dist)'), ylabel('dpos (Norm dist)')
% [r p] = corrcoef(abs(PW.doriNorm{1}),abs(PW.dposNorm{1}));
% title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])
% 
% 
% %%
% Did = 1;
% [mat xdom ydom] = smoothscatter(PW.sizesum{Did},abs(PW.dpos{Did}),.01,.01);
% 
% figure,
% subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
% xlabel('size sum'), ylabel('dpos')
% subplot(1,2,1),scatter(PW.sizesum{Did},abs(PW.dpos{Did}),'.k')
% xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
% xlabel('size sum'), ylabel('dpos')
% [r p] = corrcoef(PW.sizesum{Did},abs(PW.dpos{Did}));
% title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])
% 
% %%
% Did = 1;
% [mat xdom ydom] = smoothscatter(PW.sizesum{Did},abs(PW.dori{Did}),.01,.01);
% 
% figure,
% subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
% xlabel('size sum'), ylabel('dori')
% subplot(1,2,1),scatter(PW.sizesum{Did},abs(PW.dori{Did}),'.k')
% xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
% xlabel('size sum'), ylabel('dori')
% [r p] = corrcoef(PW.sizesum{Did},abs(PW.dori{Did}));
% title(['r = ' num2str(r(1,2))  '  p = ' num2str(p(1,2))])

function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;

