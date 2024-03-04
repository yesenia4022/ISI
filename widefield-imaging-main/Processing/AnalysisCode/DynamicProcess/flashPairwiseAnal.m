function [cc D] = flashPairwiseAnal(popResp,CoM)


cc = ones(length(CoM),length(CoM),length(popResp{1}));
D = zeros(length(CoM),length(CoM));

for i = 1:length(CoM)

    for j = i+1:length(CoM)
        
        D(i,j) = sqrt((CoM{i}(1)-CoM{j}(1))^2 + (CoM{i}(1)-CoM{j}(1))^2);

        for tau = 1:length(popResp{i})

            id = find(~isnan(popResp{i}{tau}.*popResp{j}{tau}));

            if length(id)>10
                
                R = corrcoef(popResp{i}{tau}(id),popResp{j}{tau}(id));

                cc(i,j,tau) = (R(1,2));  %cross-correlogram for cell k

                %             if cc(k) <-.6
                %                 figure,scatter(popResp{i}{tau}(id),popResp{j}{tau}(id))
                %             end


            else

                cc(i,j,tau) = NaN;
            end
            
        end

    end
end

idx = find(cc(:,:,1) ~= 1);

for i = 1:length(cc(1,1,:))
    ccdum = cc(:,:,i);
    ccdum = ccdum.*ccdum';
    cc(:,:,i) = ccdum;
end


figure
c = ceil(sqrt(length(cc(1,1,:))));
r = floor(sqrt(length(cc(1,1,:))));
for i = 1:length(cc(1,1,:))
    ccdum = cc(:,:,i);
    subplot(r,c,i)
    scatter(D(idx),ccdum(idx),'.')
    [dom mu sig] = makeEbars(D(idx),ccdum(idx),4);
    hold on
    errorbar(dom,mu,sig,'k','linewidth',4)
    
end


function [dom mu sig] = makeEbars(x,y,Nbins)

[xs id] = sort(x(:));

ys = y(id);
sampsize = floor(length(x)/Nbins);

for i = 1:Nbins
    ran = [(i-1)*sampsize + 1  i*sampsize];
    mu(i) = mean(ys(ran));
    sig(i) = std(ys(ran))/sqrt(length(ran));
    dom(i) = median(xs(ran));
end


