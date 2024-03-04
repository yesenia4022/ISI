function plotRandPosMaps

%Ian Nauhaus

global ACQinfo TC MK maskS DM cellS

[xmicperpix ymicperpix] = getImResolution;

Ncell = length(cellS.kernAll);

%Fit plane
%%
for c = 1:length(DM.colordom)


    H = [MK.CoM(:,1)*ymicperpix/1000 MK.CoM(:,2)*xmicperpix/1000 ones(Ncell,1)];

    dv = inv(H'*H)*H'*TC.ypos{c}(:);  %v slope - vertical retinotopic coordinates
    yposhat{c} = H*dv;
    du = inv(H'*H)*H'*TC.xpos{c}(:);  %u slope - horizontal retinotopic coordinates
    xposhat{c} = H*du;

    vmmperdeg = 1/sqrt(dv(1)^2 + dv(2)^2)
    ummperdeg = 1/sqrt(du(1)^2 + du(2)^2)

    dudy = du(1); dudx = du(2); dvdy = dv(1); dvdx = dv(2);  %x and y are horizontal and vertical brain coordinates
    Vratio = 1./abs(dudx*dvdy - dudy*dvdx)  %Determinant computes mm^2/deg^2

    %This way was a bigger PITA than just doing PCA
    % vecs = [dudx+1i*dudy   dvdx+1i*dvdy  -dudx-1i*dudy  -dvdx-1i*dvdy];  %each corner of the square
    % vecs2 = abs(vecs).*exp(1i*angle(vecs)*2);  %rotate
    % Res = sum(vecs2) / sum(abs(vecs))
    % distorfactor = abs(Res)  %distortion factor
    % PrincAng = angle(sqrt(Res))*180/pi  %angle of max slope (deg/mm) counter-cw from horizontal
    % prinvec = [cos(PrincAng*pi/180) sin(PrincAng*pi/180)]/sqrt(2);
    % OrthAng = angle(sqrt(Res)*exp(1i*pi/2))*180/pi; %rotate by 90
    % orthvec = [cos(OrthAng*pi/180) sin(OrthAng*pi/180)]/sqrt(2);
    % sqrt(sum((prinvec*[[dudx dudy]' [dvdx dvdy]']).^2))
    % sqrt(sum((orthvec*[[dudx dudy]' [dvdx dvdy]']).^2))

    X = [dudy dudx;  dvdy dvdx ;  -dudy -dudx ;  -dvdy -dvdx];
    [COEFF,SCORE,latent] = princomp(X);
    PrincAng = atan(COEFF(1,1)/COEFF(2,1))*180/pi
    OrthAng = atan(COEFF(1,2)/COEFF(2,2))*180/pi
    PrincMag = sqrt(latent(1))
    OrthMag = sqrt(latent(2))
    distortfactor = PrincMag/OrthMag

    dv_ang = atan(dvdy/dvdx)*180/pi;
    du_ang = atan(dudy/dudx)*180/pi;
    xyoridiff = oridiff(dv_ang*pi/180,du_ang*pi/180)*180/pi  %90 means that they are perpendicular

end



xposIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
yposIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
xposhatIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
yposhatIm = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
msk = maskS.neuronmask;
for c = 1:length(DM.colordom)    
        
    xposhat{c} = phi(xposhat{c}-prctile(TC.xpos{c},.5));
    yposhat{c} = phi(yposhat{c}-prctile(TC.ypos{c},.5));
    TC.xpos{c} = phi(TC.xpos{c}-prctile(TC.xpos{c},.5));
    TC.ypos{c} = phi(TC.ypos{c}-prctile(TC.ypos{c},.5));

    for p = 1:Ncell

        idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));
        
        xposIm(idcell) = TC.xpos{c}(p);
        yposIm(idcell) = TC.ypos{c}(p);
        
        xposhatIm(idcell) = xposhat{c}(p);
        yposhatIm(idcell) = yposhat{c}(p);

    end

    
    xyran = max([prctile(TC.xpos{c},99) prctile(TC.ypos{c},99)]);
    figure,
    subplot(2,2,1), imagesc(xposIm,'AlphaData',msk, [0 xyran]), colorbar
    title('x position'), axis image
    subplot(2,2,2), imagesc(yposIm,'AlphaData',msk, [0 xyran]), colorbar
    title('y position'), axis image
    
    subplot(2,2,3), imagesc(xposhatIm,'AlphaData',msk, [0 xyran]), colorbar
    title([num2str(ummperdeg) 'mm/deg']), axis image
    subplot(2,2,4), imagesc(yposhatIm,'AlphaData',msk, [0 xyran]), colorbar
    title([num2str(vmmperdeg) 'mm/deg']), axis image

end

    
%TC.OAng = Opreffit;
%OMag = Omagfit;

for c = 1:length(DM.colordom)

    xposRes = TC.xpos{c} - xposhat{c}';
    yposRes = TC.ypos{c} - yposhat{c}';
    doriAll = [];
    dretAll = [];
    for i = 1:length(TC.xpos{c})
        for j = i+1:length(TC.ypos{c})

            dori = oridiff(TC.OAng{c}(i)*pi/180,TC.OAng{c}(j)*pi/180)*180/pi;
            dret = sqrt((xposRes(i)-xposRes(j))^2 + (yposRes(i)-yposRes(j))^2);

            doriAll = [doriAll dori];
            dretAll = [dretAll dret];

        end
    end

    [R p] = corrcoef(doriAll,dretAll);
    R = R(1,2);
    p = p(1,2);

end
% [mat xdom ydom] = smoothscatter(doriAll,dretAll,.5,.005);
% 
% figure,
% subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
% subplot(1,2,1),scatter(doriAll,dretAll,'.k')
% xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
% title(['R = ' num2str(R)  '  p = ' num2str(p)])
% xlabel('dori (degrees)'), ylabel('dret (degrees)')


mag = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);
ang = zeros(ACQinfo.linesPerFrame,ACQinfo.pixelsPerLine);



for c = 1:length(DM.colordom)
    TC.OMag{c} = TC.OMag{c}-prctile(TC.OMag{c},0);
    TC.OMag{c} = TC.OMag{c}/prctile(TC.OMag{c},100);

    for p = 1:Ncell

        idcell = find(MK.masklabel(:) == MK.celldom(MK.nID(p)));

        mag(idcell) = TC.OMag{c}(p);
        ang(idcell) = TC.OAng{c}(p)-5;

    end

    figure,
    imagesc(ang,'AlphaData',(mag),[0 180]), colormap hsv, colorbar, axis image
    
    
    hypP = 10;
    hypO = hypP*OrthMag/PrincMag;
    hold on
    plot([0 hypP*cos(PrincAng*pi/180)]+4,[0 hypP*sin(PrincAng*pi/180)]+8,'k')
    hold on
    plot([0 hypO*cos(OrthAng*pi/180)]+4,[0 hypO*sin(OrthAng*pi/180)]+8,'k')

end