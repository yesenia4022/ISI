processF1
%Set the directories of your experiment in processF1

%Below is the approximate order of operations for getting the data saved.

% 1) load expt
% 2) Enter the repeat
% 3) process .raw 
% 4) global f1m
% 5) f1m_repX = f1m
% 6) repeat 1-5 until you have all the repeats
% 7) save('C:\saved Kmaps by repeat\xq4_002_000\','f1m_rep1','f1m_rep2','f1m_rep3')

%% Now you can load your f1 vector images (these should have been sa)
expt = 'xq4_002_000';
load(['C:\saved Kmaps by repeat\' expt])

%% filter settings
adaptbit = 0;
L = []; %no filtering
H = []; %no filtering

%%  %Make the maps from each repeat.

bw = ones(size(f1m_rep1{1}));

[kmap_hor_rep1 kmap_vert_rep1] = Gprocesskret(f1m_rep1,bw,adaptbit,L,H);
[kmap_hor_rep2 kmap_vert_rep2] = Gprocesskret(f1m_rep2,bw,adaptbit,L,H);
[kmap_hor_rep3 kmap_vert_rep3] = Gprocesskret(f1m_rep3,bw,adaptbit,L,H);

clear f1mAll
f1mAll{1} = (f1m_rep1{1} + f1m_rep2{1} + f1m_rep3{1})/3; %The mean f1 image from drift leftward, over the 3 repeats.
f1mAll{2} = (f1m_rep1{2} + f1m_rep2{2} + f1m_rep3{2})/3;
f1mAll{3} = (f1m_rep1{3} + f1m_rep2{3} + f1m_rep3{3})/3;
f1mAll{4} = (f1m_rep1{4} + f1m_rep2{4} + f1m_rep3{4})/3;
[kmap_hor_Allreps kmap_vert_Allreps] = Gprocesskret(f1mAll,bw,adaptbit,L,H);  %Generate the maps from all reapeats.


%% Make visual field sign map (see Garrett et al 2014) and use it to make your ROI

kmap_hor = kmap_hor_Allreps;  %Identify the map you want to use
kmap_vert = kmap_vert_Allreps;

kmap_hor(find(isnan(kmap_hor))) = 0;  %get rid of the NaN values.
kmap_vert(find(isnan(kmap_vert))) = 0;

sigH = 8;  %standard deviation of the Gaussian you use to smooth.
hh = fspecial('gaussian',size(kmap_hor),sigH);   %Make a gaussian filter.
hh = hh/sum(hh(:));
kmap_hor = ifft2(fft2(kmap_hor).*abs(fft2(hh)));  %Smooth hor maps with the Gaussian

hh = fspecial('gaussian',size(kmap_hor),sigH); 
hh = hh/sum(hh(:));
kmap_vert = real(ifft2(fft2(kmap_vert).*abs(fft2(hh))));  %same for vertical retinotopy map

[dhdx dhdy] = gradient(kmap_hor);  %make the gradient of the map
[dvdx dvdy] = gradient(kmap_vert);

graddir_hor = atan2(dhdy,dhdx);   %Identify the direction of the gradient for horizontal ret  (an angle value at each pixel)
graddir_vert = atan2(dvdy,dvdx);

vdiff = exp(1i*graddir_hor) .* exp(-1i*graddir_vert); %Should be vert-hor, but the gradient in Matlab for y is opposite.
VFS = sin(angle(vdiff)); %Visual field sign map (See Garrett, Nauhaus et al 2014)
id = find(isnan(VFS));
VFS(id) = 0;

hh = fspecial('gaussian',size(VFS),3); 
hh = hh/sum(hh(:));
VFS = ifft2(fft2(VFS).*abs(fft2(hh)));  %Important to smooth before thresholding below

figure,imagesc(VFS)
bw = roipoly();

%% Compute the GLOBAL coherence of the H and V retinotopy maps
kmap_hor = kmap_hor_Allreps;  %Identify the map you want to use
kmap_vert = kmap_vert_Allreps;

kmap_hor(find(isnan(kmap_hor))) = 0;  %get rid of the NaN values.
kmap_vert(find(isnan(kmap_vert))) = 0;

[dhdx dhdy] = gradient(kmap_hor);  %make the gradient of the map
[dvdx dvdy] = gradient(kmap_vert);

HorGrad = dhdx + 1i*dhdy;
VertGrad = dvdx + 1i*dvdy;

HorGrad = exp(1i*angle(HorGrad));
VertGrad = exp(1i*angle(VertGrad));

idBW = find(bw);
HorCoherenceGlobal = abs(mean(HorGrad(idBW)))
VertCoherenceGlobal = abs(mean(VertGrad(idBW)))

%% Compute the LOCAL coherence of the H and V retinotopy maps
kmap_hor = kmap_hor_Allreps;  %Identify the map you want to use
kmap_vert = kmap_vert_Allreps;

kmap_hor(find(isnan(kmap_hor))) = 0;  %get rid of the NaN values.
kmap_vert(find(isnan(kmap_vert))) = 0;

sigH = 10;  %standard deviation of the Gaussian you use to smooth.
hh = fspecial('gaussian',size(kmap_hor),sigH);   %Make a gaussian filter.
hh = hh/sum(hh(:));
normer = real(ifft2(fft2(bw).*abs(fft2(hh))));  %same for vertical retinotopy map

[dhdx dhdy] = gradient(kmap_hor);  %make the gradient of the map
[dvdx dvdy] = gradient(kmap_vert);

HorGrad = dhdx + 1i*dhdy;
VertGrad = dvdx + 1i*dvdy;

HorGrad = exp(1i*angle(HorGrad));
VertGrad = exp(1i*angle(VertGrad));

HorGrad = ifft2(fft2(HorGrad.*bw).*abs(fft2(hh)));  %Smooth hor maps with the Gaussian
%HorGrad = HorGrad./normer;

HorGrad_localCoh = abs(HorGrad);

VertGrad = real(ifft2(fft2(VertGrad.*bw).*abs(fft2(hh))));  %same for vertical retinotopy map
%VertGrad = VertGrad./normer;
VertGrad_localCoh = abs(VertGrad);

idBW = find(bw);

figure,
subplot(2,2,1), imagesc(HorGrad_localCoh,[0 .3])
hold on
title('Local coherence of horizontal retinotopy map')
subplot(2,2,2), hist(HorGrad_localCoh(idBW))
xlabel('Local coherence of horizontal retinotopy map')


subplot(2,2,3),imagesc(VertGrad_localCoh,[0 .3])
title('Local coherence of Vertical retinotopy map')
subplot(2,2,4), hist(VertGrad_localCoh(idBW))
xlabel('Local coherence of Vertical retinotopy map')

%% Compute SNR from the amplitude of the sinewave, divided be the variance of the waveform)
%% Horizontal retinotopy

SNR_amplitude_hor =  (abs(f1mAll{1})+abs(f1mAll{3}))/2;  %SNR of left and right direction
ma = prctile(SNR_amplitude_hor(:),99.9);
idma = find(SNR_amplitude_hor>=ma);
SNR_amplitude_hor(idma) = ma;
figure,imagesc(SNR_amplitude_hor,[0 ma])
title('SNR of horizontal retinotopy map')

SNR_hor_vec = SNR_amplitude_hor(find(bw(:)==1));
figure, hist(SNR_hor_vec,100)
title('SNR of horizontal retinotopy map')


%% Compute SNR from the amplitude of the sinewave, divided be the variance of the waveform)
%% Vertical retinotopy

% SNR_amplitude_vert =  (abs(f1mAll{1})+abs(f1mAll{3}))/2;  %SNR of left and right direction
% ma = prctile(SNR_amplitude_vert(:),99.9);
% idma = find(SNR_amplitude_vert>=ma);
% SNR_amplitude_vert(idma) = ma;
% figure,imagesc(SNR_amplitude_vert,[0 ma])
% title('SNR of vertical retinotopy map')
% 
% SNR_vert_vec = SNR_amplitude_vert(find(bw(:)==1));
% figure, hist(SNR_vert_vec,100)
% title('SNR of vertical retinotopy map')

