function pairWiseRCanalysis2

%Ian Nauhaus

%Generate all the pairwise values and the distance between them

%2 doesn't parse cortical distance and generate the plots

global TC MK PW DM kernC kernSigC

global sfprincax oriprincax 

[xmicperpix ymicperpix] = getImResolution;

PW = struct;

sig = 50;
dom = DM.taudom-mean(DM.taudom);
psmooth = exp(-dom.^2/(2*sig^2));
psmooth = psmooth/sum(psmooth);
psmooth = ones(length(DM.phasedom),1)*psmooth;
psmooth = abs(fft(psmooth,[],2));

Npair = 0;
for i = 1:MK.Ncell
    for j = (i+1):MK.Ncell
        Npair = Npair+1;
    end
end


for c = 1:length(TC.opref)  %length(DM.colordom) might be different 
    PW.dori{c} = zeros(1,Npair); PW.dsf{c} = zeros(1,Npair); PW.doriEuc{c} = zeros(1,Npair); PW.dsfEuc{c} = zeros(1,Npair);
    PW.Dist{c} = zeros(1,Npair); PW.ax{c} = zeros(1,Npair);
    PW.doripair{c} = zeros(Npair,2); PW.dsfpair{c} = zeros(Npair,2); 
    
    PW.dF1F0{c} = zeros(1,Npair);
    PW.dF1F0pair{c} = zeros(Npair,2);

    PW.dtcoripair{c}= zeros(Npair,2*length(TC.tcoriall_fit{c}(1,:))); PW.dtcsfpair{c} = zeros(Npair,2*length(TC.tcsfall_fit{c}(1,:)));
    PW.dtcphasepair{c} = zeros(Npair,2*length(DM.phasedom));    

    PW.dphase{c} = zeros(1,Npair);
    PW.dphasepair{c} = zeros(Npair,2);

    PW.dphaseA{c}= zeros(1,Npair);
    PW.dphaseApair{c} = zeros(Npair,2);

    k = 1;
    for i = 1:MK.Ncell

        for j = (i+1):MK.Ncell

            %Don't take abs() of dori/dsf... we want the sign to compute the
            %gradient direction
            PW.dori{c}(k) = (oridiff(TC.opref{c}(i)*pi/180,TC.opref{c}(j)*pi/180)*180/pi); %degrees
            PW.dsf{c}(k) = (log2(TC.sfpref{c}(i)/TC.sfpref{c}(j)));  %octaves
            PW.dphase{c}(k) = angle(exp(1i*pi/180*(TC.phase{c}(i)-TC.phase{c}(j))))*180/pi; %degrees
            PW.dF1F0{c}(k) = TC.F1F0{c}(i) - TC.F1F0{c}(j); 

            PW.doripair{c}(k,:) = [TC.opref{c}(i) TC.opref{c}(j)];  %useful to have the oris for later
            PW.dsfpair{c}(k,:) = [TC.sfpref{c}(i) TC.sfpref{c}(j)];  %useful to have the sf for later
            PW.dphasepair{c}(k,:) = [TC.phase{c}(i) TC.phase{c}(j)];  %useful to have the sf for later
            PW.dF1F0pair{c}(k,:) = [TC.F1F0{c}(i) TC.F1F0{c}(j)];  %useful to have the sf for later

            PW.dtcoripair{c}(k,:) = [TC.tcoriall_fit{c}(i,:) TC.tcoriall_fit{c}(j,:)];  %useful to have the oris for later
            PW.dtcsfpair{c}(k,:) = [TC.tcsfall_fit{c}(i,:) TC.tcsfall_fit{c}(j,:)];  %useful to have the oris for later

            %tcdum = TC.tcoriall{c}(i,:); tcdum2 = TC.tcoriall{c}(j,:);
            tcdum = TC.tcoriall_fit{c}(i,:); tcdum2 = TC.tcoriall_fit{c}(j,:);
            tcdum = tcdum-min(tcdum); tcdum2 = tcdum2-min(tcdum2);
            v1 = tcdum/norm(tcdum);  v2 = tcdum2/norm(tcdum2);
            cc = corrcoef(v1,v2); 
            PW.doriEuc{c}(k) = (-cc(1,2)+1)/2; %normalize between 0 and 1

            %tcdum = TC.tcsfall{c}(i,:); tcdum2 = TC.tcsfall{c}(j,:);
            tcdum = TC.tcsfall_fit{c}(i,:); tcdum2 = TC.tcsfall_fit{c}(j,:);
            tcdum = tcdum-min(tcdum); tcdum2 = tcdum2-min(tcdum2);
            v1 = tcdum/norm(tcdum);  v2 = tcdum2/norm(tcdum2);
            cc = corrcoef(v1,v2); 
            PW.dsfEuc{c}(k) = (-cc(1,2)+1)/2; %normalize between 0 and 1

            dy = (MK.CoM(i,1)-MK.CoM(j,1))*ymicperpix;
            dx = (MK.CoM(i,2)-MK.CoM(j,2))*xmicperpix;
            PW.Dist{c}(k) = sqrt(dy^2 + dx^2); %Dist between cells in microns

            PW.ax{c}(k) = atan2(dy,dx)*180/pi;  %we want atan2 to compute gradient direction
            
            if length(DM.colordom) > 1
                PW.dlumpref{c}(k) = TC.lumpref{1}(i) - TC.lumpref{1}(j);

                %PW.dLMphaseDiff{c}(k) = oridiff(TC.LMphaseDiff{1}(i)*pi/180/2,TC.LMphaseDiff{1}(j)*pi/180/2)*2*180/pi;

                PW.dLMphaseDiff{c}(k) = TC.LMphaseDiff{1}(i) - TC.LMphaseDiff{1}(j);

                PW.LMproj{c}(k) = TC.LMproj{1}(i) - TC.LMproj{1}(j);
            end

            
            %%%%Compute the pairwise for phase%%%%%%

            muori = angle(exp(1i*TC.OAng{c}(i)*pi/180) + exp(1i*TC.OAng{c}(j)*pi/180))*180/pi;
            if muori<0
                muori = muori+180;
            end

            musf = sqrt(TC.sfpref{c}(i)*TC.sfpref{c}(j));
            
            if c == 4            
                kdum_i = (kernC{1}{i} + kernC{2}{i})/2; 
                kdum_j = (kernC{1}{j} + kernC{2}{j})/2; 
                
                kSigdum_i = (kernSigC{1}{i} + kernSigC{2}{i})/2; 
                kSigdum_j = (kernSigC{1}{j} + kernSigC{2}{j})/2;
            else
                kdum_i = kernC{c}{i}; 
                kdum_j = kernC{c}{j}; 
                
                kSigdum_i = kernSigC{c}{i}; 
                kSigdum_j = kernSigC{c}{j};
            end            
            
            domdum = DM.oridom;
            if length(DM.oridom) >= 14  %downsample orientation to compute phase
                kdum_i = (kdum_i(1:2:end,:,:,:) + kdum_i(2:2:end,:,:,:))/2;
                kdum_j = (kdum_j(1:2:end,:,:,:) + kdum_j(2:2:end,:,:,:))/2;
                kSigdum_i = (kSigdum_i(1:2:end,:,:,:) + kSigdum_i(2:2:end,:,:,:))/2;
                kSigdum_j = (kSigdum_j(1:2:end,:,:,:) + kSigdum_j(2:2:end,:,:,:))/2;
                domdum = (DM.oridom(1:2:end) + DM.oridom(2:2:end))/2;
            end

            [dum oriid] = min(abs(muori-domdum));
            [dum sfid] = min(abs(musf-DM.sfdom));
            
            kernplot_i = squeeze(kdum_i(oriid,sfid,:,:));  %phase and time
            kernplot_j = squeeze(kdum_j(oriid,sfid,:,:));  %phase and time
            kernSigplot_i = squeeze(kSigdum_i(oriid,sfid,:,:));  %phase and time
            kernSigplot_j = squeeze(kSigdum_j(oriid,sfid,:,:));  %phase and time
            
            if size(kernplot_i,2) == 1  %in case there is only one phase
                kernplot_i = kernplot_i';
                kernplot_j = kernplot_j';
                
                kernSigplot_i = kernSigplot_i';
                kernSigplot_j = kernSigplot_j';                               
                
            end
            
            kernplot_i = ifft(fft(kernplot_i,[],2).*psmooth,[],2); %smooth in time
            tcphase_i = squeeze(kernplot_i(:,TC.tauID{c}(i)));
            
            kernplot_j = ifft(fft(kernplot_j,[],2).*psmooth,[],2);
            tcphase_j = squeeze(kernplot_j(:,TC.tauID{c}(j)));
            
            kernSigplot_i = ifft(fft(kernSigplot_i,[],2).*psmooth,[],2); %smooth in time
            tcphaseSig_i = squeeze(kernSigplot_i(:,TC.tauID{c}(i)));
            
            kernSigplot_j = ifft(fft(kernSigplot_j,[],2).*psmooth,[],2);
            tcphaseSig_j = squeeze(kernSigplot_j(:,TC.tauID{c}(j)));

            [ma ma_id] = max(tcphase_i); [mi mi_id] = min(tcphase_i);
            SNR_i = (ma-mi)/(tcphaseSig_i(ma_id) + tcphaseSig_i(mi_id));

            [ma ma_id] = max(tcphase_j); [mi mi_id] = min(tcphase_j);
            SNR_j = (ma-mi)/(tcphaseSig_j(ma_id) + tcphaseSig_j(mi_id));

            ma_i = mean(TC.tcphase{c}(i,:));
            ma_j = mean(TC.tcphase{c}(j,:));

            if mean(tcphase_i) > .8*ma_i & mean(tcphase_j) > .8*ma_j & TC.SNR{c}(i)>1 & TC.SNR{c}(j)>1% & SNR_i>1 & SNR_j>1

                f1_i = sum(tcphase_i'.*exp(1i*DM.phasedom*pi/180));
                f1_j = sum(tcphase_j'.*exp(1i*DM.phasedom*pi/180));

                phase_i = angle(f1_i)*180/pi;
                phase_j = angle(f1_j)*180/pi;

                PW.dphaseA{c}(k) = angle(exp(1i*(phase_i - phase_j)*pi/180))*180/pi;
                PW.dphaseApair{c}(k,:) = [phase_i phase_j];  %need the pairs for Bootstrap later
                PW.dtcphasepair{c}(k,:) = [tcphase_i' tcphase_j'];  %useful to have the oris for later

            else
                PW.dphaseA{c}(k) = NaN;
                PW.dphaseApair{c}(k,:) = [NaN NaN];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            k = k+1;

        end
    end
    
    
end

%Get color vector for each cell

%The color loop is totally redundant.  I do this only for convenience later
if length(DM.colordom) > 1
    for c = 1:length(TC.opref)
        k = 1;
        PW.dtccolorpair{c} = zeros(Npair,2*length(DM.colordom));
        PW.dcolorEuc{c} = zeros(1,Npair);
        for i = 1:MK.Ncell
            for j = (i+1):MK.Ncell

                tc1 = TC.tccolorall{c}(i,:); %Each cell 'c' is actually identical
                tc2 = TC.tccolorall{c}(j,:);

                v1 = tc1-min(tc1); v2 = tc2-min(tc2);
                v1 = v1/norm(v1);  v2 = v2/norm(v2);
                cc = corrcoef(v1,v2);
                try
                    cc = -cc(1,2);
                catch
                    'hello'
                end

                PW.dcolorEuc{c}(k) = (cc+1)/2;
                PW.dtccolorpair{c}(k,:) = [tc1 tc2];

                k = k+1;
            end
        end
    end
end

function dist = oridiff(angle1,angle2)

%pepOriDiff        Returns the difference between the angles in angle1
%                  and angle2 in the orientation domain (that is they
%                  wrap around at pi radians!!!  The angles should be in rad.

w1 = exp(1i*2*angle1);
w2 = exp(1i*2*angle2);
dist = angle(w1 ./ w2)/2;