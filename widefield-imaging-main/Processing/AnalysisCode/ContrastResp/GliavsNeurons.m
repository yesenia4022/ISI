function GliavsNeurons

global f0m maskS Analyzer

[kernPop_neuron kernPop_glia] = GetCellKernels(f0m,maskS.bwCell1,maskS.bwCell2);
kernPop{1} = kernPop_neuron;
kernPop{2} = kernPop_glia;

Nsym = length(Analyzer.loops.conds{1}.symbol);  %number of looping params
for idsym = 1:Nsym
    if strcmp('ori',Analyzer.loops.conds{1}.symbol{idsym})
        oriid = idsym;
    elseif strcmp('contrast',Analyzer.loops.conds{1}.symbol{idsym})
        contid = idsym;
    end
end

for ctype = 1:2
    dim = size(kernPop{ctype});
    
    dori = 360/dim(oriid);
    oridom = (0:dim(oriid)-1)*dori;
    
    orthD = round(90/dori)+1; %index of baseline
    slopeD1 = 2; %index for beginning of slope
    slopeD2 = 3; %index for end of slope
    
    k = 1;
    for i = 1:dim(3)   %loop through each cell
        
        kern = kernPop{ctype}(:,:,i);
        if oriid == 1
            kern = kern';
            dim(1:2) = size(kern);
        end
        
        if max(kern(:)) > .1
            
            tcdum = mean(kern(end-1:end,:));  %Take the 2 highest contrasts to get the tuning curve 'shape'
            %tcdum = kern(end,:);
            idma = find(tcdum == max(tcdum));
            kern = circshift(kern,[0 round(dim(2)/2)-idma]);
            idma = round(dim(2)/2);
            
            k1 = kern(:,idma:end);
            k2 = fliplr(kern(:,1:idma));
            L = min([length(k1(1,:)) length(k2(1,:))]);
            kern = (k1(:,1:L) + k2(:,(length(k2(1,:))-L+1):end))/2;
            
            %     id1 = find(kern(end-1,:) == max(kern(end-1,:)));
            %     id2 = find(kern(end,:) == max(kern(end,:)));
            
            %   if id1 == id2
            for c = 1:dim(1)
                base = kern(c,orthD);
                
                %kern(c,:) = kern(c,:) - base;  %Subtract response to orthogonal oris
                %
                %             kern(c,:) = kern(c,:)/kern(c,1); %Normalize the amplitude
                
                for s = 1:orthD-1
                    Mag{ctype}{c}(k,s) = log10(kern(c,s)/kern(c,s+1));  %Get slope
                end
                
            end
            
            kernPopmod{ctype}(:,:,k) = kern;
            %     else
            %         Mag{c}(k) = NaN;
            %     end
            
            k = k+1;
            
        end
    end
end

dim = size(kernPopmod{1});
oridom = oridom(1:dim(2));
% 
% mukern = mean(kernPopmod,3)';
% sigkern = std(kernPopmod,[],3)';
% ma = max(mukern);
% ma = ones(dim(2),1)*ma;
% % mukern = mukern./ma;
% % sigkern = sigkern./ma;
% mukern = mukern(:,end-1:end);
% sigkern = sigkern(:,end-1:end);
% figure,
% errorbar(oridom'*ones(1,length(mukern(1,:))),mukern,sigkern/sqrt(dim(3)));
% legend('low contrast','high contrast')
% xlabel('orientation')



% for c = 1:dim(1)    
%     id = find(Mag{end-1} > 1.2 | Mag{end-1} < 0 | Mag{end} > 1.2 | Mag{end} < 0);
%     Mag{end-1}(id) = [];
%     Mag{end}(id) = [];
% end
   



% figure,
% Nslopes = length(Mag{1}(1,:));
% for s = 1:Nslopes
%     subplot(1,Nslopes,s)
%     plot(Mag{end-1}(:,s),Mag{end}(:,s),'.'), hold on
%     plot([-1 4], [-1 4],'r')
%     xlabel('Response Ratio (low contrast)'), ylabel('Response Ratio (high contrast)')
%     hold off
%     xlim([-1 5]), ylim([-1 5])
% end



%%%%%%%%%%
%%%%%%%%%%



h = exp(1i*oridom*2*pi/180)';
for ctype = 1:2
    for i = 1:dim(1)
        kdum = squeeze(kernPopmod{ctype}(i,:,:))';
        normer = sum(kdum,2);
        pol = (kdum*h)./normer;
        Mag{ctype}{i} = abs(pol);
    end
    id = find(Mag{ctype}{end-1} > 1 | Mag{ctype}{end} > 1);
    Mag{ctype}{end-1}(id) = []; Mag{ctype}{end}(id) = [];
end

% figure,hist(Mag{end-1}-Mag{end})
% xlabel('OSI (low contrast) - OSI (high contrast)')
figure
scatter(Mag{1}{end-1},Mag{1}{end},'.'), hold on
scatter(Mag{2}{end-1},Mag{2}{end},'.r'), hold on
plot([0 1.5], [0 1.5],'k')
xlabel('OSI (low contrast)'), ylabel('OSI (high contrast)')
hold off
xlim([0 1.2]), ylim([0 1.2])
legend('neurons','glia')

figure
subplot(2,2,1)
hist(Mag{1}{end},0:.2:1), xlabel('OSI'), title('Neurons (high contrast)')
subplot(2,2,3)
hist(Mag{2}{end},0:.2:1), xlabel('OSI'), title('Glia (high contrast)')

subplot(2,2,2)
hist(Mag{1}{end-1},0:.2:1), xlabel('OSI'), title('Neurons (low contrast)')
subplot(2,2,4)
hist(Mag{2}{end-1},0:.2:1), xlabel('OSI'), title('Glia (low contrast)')




    