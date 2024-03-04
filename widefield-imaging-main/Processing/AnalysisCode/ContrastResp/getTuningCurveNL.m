function [dum1 dum2] = getTuningCurveNL

%kernPop can be obtained from 'GetCellKernels'

global f0m maskS

kernPop = GetCellKernels(f0m,maskS.bwCell1,maskS.bwCell2);

dim = size(kernPop);

dori = 360/dim(2);
oridom = (0:dim(2)-1)*dori;

orthD = round(90/dori)+1; %index of baseline

k = 1;
for i = 1:dim(3)   %loop through each cell

    kern = kernPop(:,:,i);

    if max(kern(:)) > .0

        tcdum = mean(kern(end-1:end,:));  %Take the 2 highest contrasts to get the tuning curve 'shape'
        %tcdum = kern(end,:);
        idma = find(tcdum == max(tcdum));
        kern = circshift(kern,[0 round(dim(2)/2)-idma]);
        idma = round(dim(2)/2);

        for c = 1:dim(1)
            base = (kern(c,idma-orthD+1) + kern(c,idma+orthD-1))/2;

            %kern(c,:) = kern(c,:) - base;  %Subtract response to orthogonal oris
            %kern(c,:) = kern(c,:)/kern(c,idma); %Normalize the amplitude
            %kern(c,:) = kern(c,:)/norm(kern(c,:));

        end
        

        %kernPopmod(:,:,k) = kern(:,idma-orthD+1:idma+orthD-1);
        kernPopmod(:,:,k) = (kern/norm(kern(:)));
        %kernPopmod(:,:,k) = kern;

        k = k+1;

    end

end

dum1 = squeeze(kernPopmod(2,:,:));
dum2 = squeeze(kernPopmod(3,:,:));






% [param] = FtoSpikePoly(dum1(:),dum2(:));
% 
% figure,plot(dum1(:),dum2(:),'.')
% 
% 
% %%%
% %%%
% beta = param(1);
% A = param(2);
% gam = param(3);
% K = param(4);
% 
% 
% 
% 
% 
% F1 = 0:.001:K;
% 
% 
% %F1 = K./(1+exp(-beta*A*gam - A*log(K./F2-1)-gam));
% 
% s1 = gam - 1./beta * log(K./F1 - 1);
% s2 = A*s1;
% F2 = K./(1+exp(-beta*(s2-gam)));
% 
% figure,
% subplot(2,2,3)
% plot(dum1(:),dum2(:),'.'), xlim([0 K]), ylim([0 K])
% hold on
% plot(F1,F2,'k', 'lineWidth',2)
% xlabel('Ca2++ Fluorescence (low contrast)'), ylabel('Ca2++ Fluorescence (high contrast)')
% subplot(2,2,4)
% plot(s1,F1,'k', 'lineWidth',2)
% xlabel('SpikeRate (high contrast)'), xlim([min(s1) max(s1)])
% subplot(2,2,1)
% plot(F2,s2,'k', 'lineWidth',2)
% ylabel('SpikeRate (low contrast)'), ylim([min(s2) max(s2)])
% subplot(2,2,2)
% plot(s1,s2,'k', 'lineWidth',2),ylim([min(s2) max(s2)]), xlim([min(s1) max(s1)])
