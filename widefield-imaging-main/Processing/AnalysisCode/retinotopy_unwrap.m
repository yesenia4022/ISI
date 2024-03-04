function [Hor_mod Vert_mod] = retinotopy_unwrap(H_pref_sphase,V_pref_sphase,xyLocations)


%load '/Users/in2293/Desktop/Random Matlab code/Data/xt7/xt7_004_001/fits.mat'
%load '/Users/in2293/Desktop/Random Matlab code/Data/xt7/xt7_004_001/xyLocations.mat'

%load '/Users/in2293/Desktop/Random Matlab code/Data/xx3/xx3_001_005/fits.mat'
%load '/Users/in2293/Desktop/Random Matlab code/Data/xx3/xx3_001_005/xyLocations.mat'

%[Hor_mod Vert_mod] = retinotopy_unwrap(fits.Hor.pref_sphase,fits.Vert.pref_sphase,xyLocations)


%%

Hor_vec = exp(H_pref_sphase(1:end)*1i*pi/180); %convert phase to vector
Vert_vec = exp(V_pref_sphase(1:end)*1i*pi/180);

mu_hor = mean(Hor_vec);  %Get the resultant
mu_vert = mean(Hor_vec);

Hor_vec = Hor_vec*conj(mu_hor)*exp(1i*pi);  %subtract the resultant and add pi
Vert_vec = Vert_vec*conj(mu_vert)*exp(1i*pi);

Hor_mod = angle(Hor_vec)*180/pi;   %convert back to phase
Vert_mod = angle(Vert_vec)*180/pi;

%Kmeans  (This did not work as well)
% opts = statset('Display','final');
% idx_H = kmeans(zscore([Hor_mod(:) xyLocations(:,:)]),2,'Distance','cityblock',...
%     'Replicates',5,'Options',opts);
% idx_V = kmeans(zscore([Vert_mod(:) xyLocations(:,:)]),2,'Distance','cityblock',...
%     'Replicates',5,'Options',opts);


%Gaussian mixture model (Same goal but better than K-means)
C = fitgmdist(zscore([Vert_mod(:) xyLocations(:,:)]),2,'Replicates',50,'RegularizationValue',0.001);
idx_V = cluster(C,zscore([Vert_mod(:) xyLocations(:,:)]));
C = fitgmdist(zscore([Hor_mod(:) xyLocations(:,:)]),2,'Replicates',50,'RegularizationValue',0.001);
idx_H = cluster(C,zscore([Hor_mod(:) xyLocations(:,:)]));

idH1 = find(idx_H == 1);  %ID of each cluster
idH2 = find(idx_H == 2);
idV1 = find(idx_V == 1);
idV2 = find(idx_V == 2);



figure,
subplot(4,2,1)
scatter(xyLocations(idH1,1),Hor_mod(idH1),'.'), hold on
scatter(xyLocations(idH2,1),Hor_mod(idH2),'.'), 
ylabel('horizontal')
subplot(4,2,3)
scatter(xyLocations(idH1,2),Hor_mod(idH1),'.'), hold on
scatter(xyLocations(idH2,2),Hor_mod(idH2),'.'), 
ylabel('horizontal')
subplot(4,2,5)
scatter(xyLocations(idV1,1),Vert_mod(idV1),'.'), hold on
scatter(xyLocations(idV2,1),Vert_mod(idV2),'.'), 
ylabel('vertical')
subplot(4,2,7)
scatter(xyLocations(idV1,2),Vert_mod(idV1),'.'), hold on
scatter(xyLocations(idV2,2),Vert_mod(idV2),'.'), 
ylabel('vertical')


%Unwrap horizontal: find cluster with smaller values, on average, and add 360
mu1 = mean(Hor_mod(idH1));
mu2 = mean(Hor_mod(idH2));
if mu1<mu2
    Hor_mod(idH1) = Hor_mod(idH1)+360;
else
    Hor_mod(idH2) = Hor_mod(idH2)+360;
end

%Unwrap horizontal: find cluster with smaller values, on average, and add 360
mu1 = mean(Vert_mod(idV1));
mu2 = mean(Vert_mod(idV2));
if mu1<mu2
    Vert_mod(idV1) = Vert_mod(idV1)+360;
else
    Vert_mod(idV2) = Vert_mod(idV2)+360;
end

subplot(4,2,2)
scatter(xyLocations(:,1),Hor_mod,'.'), 
ylabel('horizontal unwrapped')
subplot(4,2,4)
scatter(xyLocations(:,2),Hor_mod,'.'), 
ylabel('horizontal unwrapped')
subplot(4,2,6)
scatter(xyLocations(:,1),Vert_mod,'.'), 
ylabel('vertical unwrapped')
subplot(4,2,8)
scatter(xyLocations(:,2),Vert_mod,'.'),
ylabel('vertical unwrapped')
