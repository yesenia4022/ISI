function localMapSlope

%Ian Nauhaus

%% Local neighborhood analysis

global TC DM MK NB

NB = struct;  %make sure to reset

[xmicperpix ymicperpix] = getImResolution;

R = 30;
for c = 1:length(DM.colordom)
    %First normalize each tuning curve
    tcoriallmod = max(TC.tcoriall{c}')'*ones(1,length(TC.tcoriall{c}(1,:)));
    tcoriall = TC.tcoriall{c}./tcoriallmod;
    tcsfallmod = max(TC.tcsfall{c}')'*ones(1,length(TC.tcsfall{c}(1,:)));
    tcsfall = TC.tcsfall{c}./tcsfallmod;
    

    for p = 1:MK.Ncell

        dy = (MK.CoM(p,1) - MK.CoM(:,1))*ymicperpix;
        dx = (MK.CoM(p,2) - MK.CoM(:,2))*xmicperpix;
        r = sqrt(dx.^2+dy.^2); %all distances in microns        
               
        idR = find(r<R);
        
        H = [MK.CoM(idR,:) ones(length(idR),1)];
        sfp = TC.sfpref{c}(idR); 
        orip = TC.opref{c}(idR);

        id = find(~isnan(log2(sfp)) & ~isinf(log2(sfp)));
        sfplane = inv(H(id,:)'*H(id,:))*H(id,:)'*log2(sfp(id))';
        sfpref_hat = 2.^(H*sfplane);
        
        shat = log(sfpref_hat(id)');
        s = log(sfp(id));
        %varacc = (var(s)-var(s-shat))/var(s);        
        
        if length(s)>5
            rsf = corrcoef(s,shat);
            if ~isnan(rsf(1,1))
                rsf = rsf(1,2);
            end
        else
            rsf = NaN;
        end

        id = find(~isnan(orip));
        H = H(id,:); orip = orip(id);
        [oriplane oripref_hat] = Planefit(H(:,2),H(:,1),orip);
        
        rori = circCorr(orip(:)*pi/180,oripref_hat(:)*pi/180);
        
        if rori>.5 & rsf>.4

            ori_ang = atan(oriplane(2)/oriplane(1))*180/pi;
            sf_ang = atan(sfplane(1)/sfplane(2))*180/pi;
            NB.maporidiff{c}(p) = abs(oridiff(ori_ang*pi/180,sf_ang*pi/180)*180/pi);

            mapslopemag{c}(p) = sqrt(oriplane(1)^2 + oriplane(2)^2) * sqrt(sfplane(1)^2 + sfplane(2)^2);

        else
            
             mapslopemag{c}(p) = NaN;
             NB.maporidiff{c}(p) = NaN;

        end
        

        
    end
    figure,hist(NB.maporidiff{c},[0:10:90])
end



