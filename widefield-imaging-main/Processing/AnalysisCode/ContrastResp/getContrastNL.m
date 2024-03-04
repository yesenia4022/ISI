function getContrastNL

%kernPop can be obtained from 'GetCellKernels'

global f0m maskS

kernPop = GetCellKernels(f0m,maskS.bwCell1,maskS.bwCell2);


dim = size(kernPop);

kernPopmod = zeros(dim);
k = 1;
for i = 1:dim(3)   %loop through each cell    
    
    kern = kernPop(:,:,i);
    
    if max(kern(:)) > .05

        tcdum = mean(kern(2:end,:));
        %tcdum = kern(end,:);
        idma = find(tcdum == max(tcdum));
        kern = circshift(kern,[0 1-idma]);

        kernPopmod(:,:,k) = kern;


        for c = 1:dim(1)

            Slope(c,k) = log((kern(c,1)/kern(c,2)));  %Get slope

        end
        
        k = k+1;
    end

    
end


normer = ones(length(Slope(:,1)),1)*max(Slope(1:end,:)); 
Slope = Slope./normer;

id = find(Slope > 10 | Slope < -1);
Slope(id) = NaN;

muSlope = nanmean(Slope,2);
sigSlope = nanstd(Slope,[],2)/sqrt(length(Slope(1,:)));

figure,errorbar(muSlope,sigSlope,'-o')

figure,plot(Slope)


