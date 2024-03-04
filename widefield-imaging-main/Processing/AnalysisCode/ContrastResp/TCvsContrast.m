function TCvsContrast

global f0m maskS Analyzer

[kernPop kernPop_glia] = GetCellKernels(f0m,maskS.bwCell1,maskS.bwCell2);
d = size(kernPop); dg = size(kernPop_glia);
kernPopdum = zeros(d(1),d(2),d(3)+dg(3));
kernPopdum(:,:,1:d(3)) = kernPop;
kernPopdum(:,:,d(3)+1:end) = kernPop_glia;
kernPop = kernPopdum;

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping params
for idsym = 1:Nsym
    if strcmp('ori',Analyzer.loops.conds{1}.symbol{idsym})
        oriid = idsym;
    elseif strcmp('contrast',Analyzer.loops.conds{1}.symbol{idsym})
        contid = idsym;
    end
end

dim = size(kernPop);

oridom = eval(Analyzer.L.param{oriid}{2});
dori = oridom(2)-oridom(1);

orthD = round(90/dori)+1; %index of baseline
slopeD1 = 2; %index for beginning of slope
slopeD2 = 3; %index for end of slope

smoother = zeros(1,length(oridom));
smoother(1:3) = [.2 1 .2];
smoother = smoother/sum(smoother);
smoother = abs(fft(smoother));

k = 1;
for i = 1:dim(3)   %loop through each cell

    kern = kernPop(:,:,i);
%     id = find(kern<0);
%     kern(id) = 0;


    if oriid == 1
        kern = kern';
        dim(1:2) = size(kern);
    end

    if max(kern(2,:)) > .05
        %if abs(idhi-idlo)<=0

        for c = 1:dim(1)  %each contrast

            kernC = squeeze(kern(c,:));           

            %kernC = ifft(fft(kernC).*smoother);

            tcdum = kernC;
            [dum idma] = max(tcdum);
            kernC = circshift(kernC,[0 1-idma]);

            k1 = kernC(1:end/2+1);
            k2 = [kernC(end/2+1:end) kernC(1)];
            kernC_peak = [kernC(end-orthD+1:end) kernC(1:orthD)];
            kernC_peak2 = kernC(orthD-2:end-orthD+1);
            %kernC_peak = kernC_peak(2:end-1) + kernC_peak2(2:end-1);

            k2 = fliplr(k2);
            kernC_wrap = (k1 + k2)/2;
            
            %Use max contrast tuning curve to determine orientation around
            %which gain values are computed

            [BW{c}(k) tcI domI idI] = FWHM(kernC_wrap(1:orthD),dori);

%            oridomI = linspace(0,oridom(length(kernC_peak)),100);            
%             kernC_peakI = interp1(oridom(1:length(kernC_peak)),kernC_peak,oridomI);
%             [param ffit varacc] = Gaussfit(1,kernC_peakI,1);
%             doriI = oridomI(2)-oridomI(1);
% % 
%             if varacc > .5 && param(2)*doriI < 80
%                 BW{c}(k) = param(2)*doriI;
%                 %BW{c}(k) = param(3);
%                 %figure,plot(kernC_peak)
%             else
%                 BW{c}(k) = NaN;
%             end

            Mag{c}(k,1) = log10(kernC_wrap(1)/kernC_wrap(2));
            Mag{c}(k,2) = log10(kernC_wrap(2)/kernC_wrap(4));

%             Mag{c}(k,1) = (kernC_wrap(1)-kernC_wrap(2))/(kernC_wrap(1)+kernC_wrap(2));
%             Mag{c}(k,2) = (kernC_wrap(2)-kernC_wrap(4))/(kernC_wrap(2))+kernC_wrap(4);

            %            base = kern(c,orthD);

            %kern(c,:) = kern(c,:) - base;  %Subtract response to orthogonal oris
            %
            %             kern(c,:) = kern(c,:)/kern(c,1); %Normalize the amplitude

            %             for s = 1:orthD-1
            %                 Mag{c}(k,s) = log10(kern(c,s)/kern(c,s+1));  %Get slope
            %             end

            kernPopmod_full(c,:,k) = kern(c,:);
            kernPopmod(c,:,k) = kernC_wrap;
            kernPopmod_pk(c,:,k) = abs(diff(kernC_peak));
            
        end


        %     else
        %         Mag{c}(k) = NaN;
        %     end

        k = k+1;

    end

end

dum1 = kernPopmod_pk(1,:,:);
dum2 = kernPopmod_pk(2,:,:);

figure,scatter(dum1(:),dum2(:),'.')
hold on
plot([0 .5],[0 .5])
xlabel('low contrast')
ylabel('high contrast')

figure,plot(BW{1},BW{2},'.'), hold on, plot([0 150],[0 150],'k')
xlabel('FWHM (low contrast)'), ylabel('FWHM (high contrast)')

dim = size(kernPopmod);
oridomTrunc = oridom(1:dim(2));

mukern = mean(kernPopmod,3)';
sigkern = std(kernPopmod,[],3)';
ma = max(mukern);
ma = ones(dim(2),1)*ma;
% mukern = mukern./ma;
% sigkern = sigkern./ma;
mukern = mukern(:,end-1:end);
sigkern = sigkern(:,end-1:end);
figure,
errorbar(oridomTrunc'*ones(1,length(mukern(1,:))),mukern,sigkern/sqrt(dim(3)));
legend('low contrast','high contrast')
xlabel('orientation')



% for c = 1:dim(1)
%     id = find(Mag{end-1} > 1.2 | Mag{end-1} < 0 | Mag{end} > 1.2 | Mag{end} < 0);
%     Mag{end-1}(id) = [];
%     Mag{end}(id) = [];
% end




figure,
Nslopes = length(Mag{1}(1,:));
for s = 1:Nslopes
    subplot(1,Nslopes,s)
    plot(Mag{end-1}(:,s),Mag{end}(:,s),'.'), hold on
    plot([-1 4], [-1 4],'r')
    xlabel('Response Ratio (low contrast)'), ylabel('Response Ratio (high contrast)')
    hold off
    xlim([-1 5]), ylim([-1 5])
end



%%%%%%%%%%
%%%%%%%%%%

% 
h = exp(1i*oridom*2*pi/180)';
for i = 1:dim(1)
    kdum = squeeze(kernPopmod_full(i,:,:))';
    normer = sum((kdum),2);
    pol = (kdum*h)./normer;
    Mag{i} = abs(pol);
end
id = find(Mag{end-1} > 1 | Mag{end} > 1);
Mag{end-1}(id) = []; Mag{end}(id) = [];

% clear Mag
% p = 1;
% h = exp(1i*oridom*2*pi/180)';
% for i = 1:length(kernPop(1,1,:))  
%     
%     kern = kernPop(:,:,i);
%     if max(kern(:))>.0
%         for k = 1:2
%             
%             tc = squeeze(kern(:,end-k+1))';
%             
%             pol = (tc*h)./sum(tc);
%             Mag(p,k) = abs(pol);
%         end
%         p = p+1;
%     end
% end
% id = find(Mag{end-1} > 1 | Mag{end} > 1);
% Mag{end-1}(id) = []; Mag{end}(id) = [];


figure,scatter(Mag{end-1},Mag{end},'.'), hold on
plot([0 1.5], [0 1.5],'r')
xlabel('OSI (low contrast)'), ylabel('OSI (high contrast)')
hold off
xlim([0 1.2]), ylim([0 1.2])

%figure,hist(Mag{end})


function [BW tcI domI id] = FWHM(tc,dori)

I = 20;

dom = linspace(0,dori*(length(tc)-1),length(tc));
domI = linspace(0,dori*(length(tc)-1),length(tc)*I);

tcI = interp1(dom,tc,domI);

%tcdum = tcI-min(tc);
tcdum = tcI;
tcdum = tcdum/max(tcdum);

%tcdum = tcI/max(tcI)

[dum id] = min(abs(tcdum-1/sqrt(2)));

BW = 2*domI(id)+ randn(1);




    