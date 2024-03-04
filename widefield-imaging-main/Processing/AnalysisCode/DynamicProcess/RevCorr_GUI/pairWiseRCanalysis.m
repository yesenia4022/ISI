function pairWiseRCanalysis

%Ian Nauhaus

global TC MK PW DM kernC kernSigC

global sfprincax oriprincax 

[xmicperpix ymicperpix] = getImResolution;

sig = 50;
dom = DM.taudom-mean(DM.taudom);
psmooth = exp(-dom.^2/(2*sig^2));
psmooth = ones(4,1)*psmooth;
psmooth = abs(fft(psmooth,[],2));

Npair = 0;
for i = 1:MK.Ncell
    for j = (i+1):MK.Ncell
        Npair = Npair+1;
    end
end

doriAll = zeros(1,Npair); dsfAll = zeros(1,Npair); doriEucAll = zeros(1,Npair); dsfEucAll = zeros(1,Npair); 
DistAll = zeros(1,Npair); axAll = zeros(1,Npair); 
doripairAll = zeros(Npair,2); dsfpairAll = zeros(Npair,2);

dtcoripairAll = zeros(Npair,2*length(TC.tcoriall_fit{1}(1,:))); dtcsfpairAll = zeros(Npair,2*length(TC.tcsfall_fit{1}(1,:)));
dtcphasepairAll = zeros(Npair,2*4);

dphaseAll = zeros(1,Npair);
dphasepairAll = zeros(Npair,2);

dphaseAAll = zeros(1,Npair);
dphaseApairAll = zeros(Npair,2);

k = 1;
for i = 1:MK.Ncell
    
    for j = (i+1):MK.Ncell
        
        %Don't take abs() of dori/dsf... we want the sign to compute the
        %gradient direction
        dori = (oridiff(TC.opref{1}(i)*pi/180,TC.opref{1}(j)*pi/180)*180/pi); %degrees
        dsf = (log2(TC.sfpref{1}(i)/TC.sfpref{1}(j)));  %octaves
        dphase = angle(exp(1i*pi/180*(TC.phase{1}(i)-TC.phase{1}(j))))*180/pi; %degrees
        doriAll(k) = dori;
        dsfAll(k) = dsf;
        dphaseAll(k) = dphase;
        doripairAll(k,:) = [TC.opref{1}(i) TC.opref{1}(j)];  %useful to have the oris for later
        dsfpairAll(k,:) = [TC.sfpref{1}(i) TC.sfpref{1}(j)];  %useful to have the sf for later
        dphasepairAll(k,:) = [TC.phase{1}(i) TC.phase{1}(j)];  %useful to have the sf for later
        
        dtcoripairAll(k,:) = [TC.tcoriall_fit{1}(i,:) TC.tcoriall_fit{1}(j,:)];  %useful to have the oris for later
        dtcsfpairAll(k,:) = [TC.tcsfall_fit{1}(i,:) TC.tcsfall_fit{1}(j,:)];  %useful to have the oris for later
        
        %tcdum = TC.tcoriall{1}(i,:); tcdum2 = TC.tcoriall{1}(j,:);
        tcdum = TC.tcoriall_fit{1}(i,:); tcdum2 = TC.tcoriall_fit{1}(j,:);
        tcdum = tcdum-min(tcdum); tcdum2 = tcdum2-min(tcdum2);
        v1 = tcdum/norm(tcdum);  v2 = tcdum2/norm(tcdum2);
        %doriEuc = 1-v1(:)'*v2(:);
        %doriEuc = norm(v1-v2)/2;  
        doriEuc = corrcoef(v1,v2); doriEuc = -doriEuc(1,2); doriEuc = (doriEuc+1)/2;
        
        %tcdum = TC.tcsfall{1}(i,:); tcdum2 = TC.tcsfall{1}(j,:);
        tcdum = TC.tcsfall_fit{1}(i,:); tcdum2 = TC.tcsfall_fit{1}(j,:);
        tcdum = tcdum-min(tcdum); tcdum2 = tcdum2-min(tcdum2);
        v1 = tcdum/norm(tcdum);  v2 = tcdum2/norm(tcdum2);        
        %dsfEuc = 1-v1(:)'*v2(:);
        %dsfEuc = norm(v1-v2)/2; 
        dsfEuc = corrcoef(v1,v2); dsfEuc = -dsfEuc(1,2); dsfEuc = (dsfEuc+1)/2;
        
        doriEucAll(k) = doriEuc;
        dsfEucAll(k) = dsfEuc;
        
        dy = (MK.CoM(i,1)-MK.CoM(j,1))*ymicperpix;
        dx = (MK.CoM(i,2)-MK.CoM(j,2))*xmicperpix;
        Dist = sqrt(dy^2 + dx^2); %Dist between cells in microns
        DistAll(k) = Dist;
        
        axAll(k) = atan2(dy,dx)*180/pi;  %we want atan2 to compute gradient direction        
        
        
        %%%%Compute the pairwise for phase%%%%%%

        muori = angle(exp(1i*TC.OAng{1}(i)*pi/180) + exp(1i*TC.OAng{1}(j)*pi/180))*180/pi;
        if muori<0
            muori = muori+180;
        end

        musf = sqrt(TC.sfpref{1}(i)*TC.sfpref{1}(j));

        [dum oriid] = min(abs(muori-DM.oridom));
        [dum sfid] = min(abs(musf-DM.sfdom));

        kernplot_i = squeeze(kernC{1}{i}(oriid,sfid,:,:));  %phase and time
        kernplot_i = ifft(fft(kernplot_i,[],2).*psmooth,[],2); %smooth in time
        tcphase_i = squeeze(kernplot_i(:,TC.tauID{1}(i)));

        kernplot_j = squeeze(kernC{1}{j}(oriid,sfid,:,:));  %phase and time
        kernplot_j = ifft(fft(kernplot_j,[],2).*psmooth,[],2);
        tcphase_j = squeeze(kernplot_j(:,TC.tauID{1}(j)));
        
        kernSigplot_i = squeeze(kernSigC{1}{i}(oriid,sfid,:,:));  %phase and time
        kernSigplot_i = ifft(fft(kernSigplot_i,[],2).*psmooth,[],2); %smooth in time
        tcphaseSig_i = squeeze(kernSigplot_i(:,TC.tauID{1}(i)));

        kernSigplot_j = squeeze(kernSigC{1}{j}(oriid,sfid,:,:));  %phase and time
        kernSigplot_j = ifft(fft(kernSigplot_j,[],2).*psmooth,[],2);
        tcphaseSig_j = squeeze(kernSigplot_j(:,TC.tauID{1}(j)));
        
        [ma ma_id] = max(tcphase_i); [mi mi_id] = min(tcphase_i);
        SNR_i = (ma-mi)/(tcphaseSig_i(ma_id) + tcphaseSig_i(mi_id));
        
        [ma ma_id] = max(tcphase_j); [mi mi_id] = min(tcphase_j);
        SNR_j = (ma-mi)/(tcphaseSig_j(ma_id) + tcphaseSig_j(mi_id));
      
        ma_i = mean(TC.tcphase{1}(i,:));
        ma_j = mean(TC.tcphase{1}(j,:));

        if mean(tcphase_i) > .8*ma_i & mean(tcphase_j) > .8*ma_j & TC.SNR{1}(i)>1 & TC.SNR{1}(j)>1 & SNR_i>1 & SNR_j>1

            f1_i = sum(tcphase_i'.*exp(1i*DM.phasedom*pi/180));
            f1_j = sum(tcphase_j'.*exp(1i*DM.phasedom*pi/180));

            phase_i = angle(f1_i)*180/pi;
            phase_j = angle(f1_j)*180/pi;

            dphaseAAll(k) = angle(exp(1i*(phase_i - phase_j)*pi/180))*180/pi;            
            dphaseApairAll(k,:) = [phase_i phase_j];  %need the pairs for Bootstrap later
            dtcphasepairAll(k,:) = [tcphase_i' tcphase_j'];  %useful to have the oris for later
            
        else
            dphaseAAll(k) = NaN;
            dphaseApairAll(k,:) = [NaN NaN]; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        k = k+1;
        
    end
end


Ddom = 0:100:300;
%Ddom = [0 40 92 160 252]; %logarithmic spacing
clear dori dsf dphase dphaseA doriEuc dsfEuc
for i = 1:length(Ddom)-1

    DLimitMin = Ddom(i);  %Limit pairs to be this distance (microns) apart
    DLimitMax = Ddom(i+1);  %Limit pairs to be this distance (microns) apart
    id = find(DistAll<DLimitMax & DistAll>DLimitMin & ~isnan(doriAll.*dsfAll) & ~isinf(doriAll.*dsfAll));

    dori{i} = doriAll(id);
    dsf{i} = dsfAll(id);
    dphase{i} = dphaseAll(id);
    dphaseA{i} = dphaseAAll(id);

    doriEuc{i} = doriEucAll(id);
    dsfEuc{i} = dsfEucAll(id);
    
    ax{i} = axAll(id);
    
    dist{i} = DistAll(id);
    
    doripair{i} = doripairAll(id,:);
    dsfpair{i} = dsfpairAll(id,:);
    dphasepair{i} = dphasepairAll(id,:);
    dphaseApair{i} = dphaseApairAll(id,:);
    
    dtcoripair{i} = dtcoripairAll(id,:);
    dtcsfpair{i} = dtcsfpairAll(id,:);
    dtcphasepair{i} = dtcphasepairAll(id,:);

end

%% Get preferred map axis
varflag = 1;  %if we want it to be non-directional (axis), set to 1
f = 1+varflag;

%I already do this for all ROIs in sfOriaxes.  I just put it here so that I
%could see it with the maps
axAll = [axAll axAll+180];
axAll = angle(exp(1i*axAll*pi/180*f))*180/pi/f;
doriAll = [doriAll -doriAll];
dsfAll = [dsfAll -dsfAll];
DistAll = [DistAll DistAll];

if varflag
    doriAll = abs(doriAll);
    dsfAll = abs(dsfAll);
end

id = find(DistAll<150);
DistAll = DistAll(id); axAll = axAll(id); doriAll = doriAll(id); dsfAll = dsfAll(id);

axmax = 90*(2-varflag);
axedges = -axmax:10:axmax; dax = axedges(2)-axedges(1);
axW = 2*dax;
axdom = axedges(1:end-1)+dax/2;
doriAlln = doriAll./DistAll*1000; %deg/mm
dsfAlln = dsfAll./DistAll*1000; %oct/mm


for j = 1:length(axdom)  %loop through axis domain

    id = find(axAll > axdom(j)-axW/2 & axAll < axdom(j)+axW/2 );

    dori_mu(j) = trimmean(doriAlln(id),20);
    dori_SE(j) = std(doriAlln(id))/sqrt(length(id));

    dsf_mu(j) = trimmean(dsfAlln(id),20);
    dsf_SE(j) = std(dsfAlln(id))/sqrt(length(id));

end

% [mat xdom ydom] = smoothscatter(axAll,doriAlln,8,20,[-180 180],[prctile(doriAlln,2) prctile(doriAlln,98)]);
% [dum id] = max(mat);
% dori_mu = ydom(id);
% 
% id = find(dsfAlln~=0);
% [mat xdom ydom] = smoothscatter(axAll(id),dsfAlln(id),4,1,[-180 180],[prctile(dsfAlln,2) prctile(dsfAlln,98)]);
% [dum id] = max(mat);
% dsf_mu = ydom(id);
% axdom = xdom;


Resori = sum(dori_mu.*exp(1i*axdom*pi/180*f))/(.5*length(dori_mu));
oriprincax = angle(Resori)*180/pi/f;
Ressf = sum(dsf_mu.*exp(1i*axdom*pi/180*f))/(.5*length(dsf_mu));
sfprincax = angle(Ressf)*180/pi/f;
M = abs(Ressf.*Resori);
daxis = abs(oridiff(sfprincax*pi/180,oriprincax*pi/180)*180/pi);

figure,
subplot(2,1,1)
scatter(axAll,doriAlln,'.')
hold on
%errorbar(axdom,dori_mu,dori_SE,'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
plot(axdom,dori_mu,'k','LineWidth',3), ylim([prctile(doriAlln,2) prctile(doriAlln,98)])
title(['phase = ' num2str(oriprincax) ' amp = ' num2str(abs(Resori))])
ylabel('deg/mm')

subplot(2,1,2)
scatter(axAll,dsfAlln,'.')
hold on
%errorbar(axdom,dsf_mu,dsf_SE,'k','LineWidth',3), ylim([prctile(dsfAlln,2) prctile(dsfAlln,98)])
plot(axdom,dsf_mu,'k','LineWidth',3), ylim([prctile(dsfAlln,2) prctile(dsfAlln,98)])
title(num2str(sfprincax))
title(['phase = ' num2str(sfprincax) ' amp = ' num2str(abs(Ressf)) '  diff = ' num2str(daxis)])
ylabel('oct/mm')
%%
[mat xdom ydom] = smoothscatter(abs(dori{1}),abs(dsf{1}),.8,.05);


figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori'), ylabel('dsf')
subplot(1,2,1),scatter(abs(dori{1}),abs(dsf{1}),'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori'), ylabel('dsf')

%%
[mat xdom ydom] = smoothscatter(doriEuc{1},dsfEuc{1},.015,.015);

figure,
subplot(1,2,2),imagesc(xdom,ydom,mat), axis xy
xlabel('dori (Euc dist)'), ylabel('dsf (Euc dist)')
subplot(1,2,1),scatter(doriEuc{1},dsfEuc{1},'.k')
xlim([xdom(1) xdom(end)]), ylim([ydom(1) ydom(end)])
xlabel('dori (Euc dist)'), ylabel('dsf (Euc dist)')

PW.dori = dori;
PW.dsf = dsf;
PW.dphase = dphase;
PW.dphaseA = dphaseA;
PW.Ddom = Ddom;
PW.doriEuc = doriEuc;
PW.dsfEuc = dsfEuc;
PW.ax = ax;
PW.dist = dist;
PW.doripair = doripair;
PW.dsfpair = dsfpair;
PW.dphasepair = dphasepair;
PW.dphaseApair = dphaseApair;
PW.dtcoripair = dtcoripair;
PW.dtcsfpair = dtcsfpair;
PW.dtcphasepair = dtcphasepair;

function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2