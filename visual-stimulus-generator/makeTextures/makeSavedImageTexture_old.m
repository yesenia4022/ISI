function makeSavedImageTexture

%This one builds only the images that get played for this trial.
%Importantly, this one doesn't allow for drift or contrast reversal like
%the old version in the 'archive' file.

global Mstate screenNum loopTrial %movieBlock

global Gtxtr TDim  %'playgrating' will use these

Screen('Close')  %First clean up: Get rid of all textures/offscreen windows

Gtxtr = []; TDim = [];  %reset

P = getParamStruct;

screenRes = Screen('Resolution',screenNum); %sometimes this reports 0Hz
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

%Replace it with this (comes from flipInterval)
screenRes.hz = Mstate.refresh_rate;



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Create Sequence for this trial%%%%%

%colordom = getColorDomain(P.colorspace);

N_Im = round(P.stim_time*screenRes.hz/P.h_per); %number of images to present


%%
stimulusFolder = P.image_folder; %'/Users/iannauhaus/Desktop/Textures';

trialList = dir(fullfile(stimulusFolder, P.fileID)); %Conditioned list in folder

imdom = trialList;  %image "domain" is the list of potential image file names


%%

immesh = ndgrid(1:length(imdom)); %degenerate use of ndgrid, but keeping consistent with other modules
%[phasemesh orimesh sfmesh colormesh] = ndgrid(1:length(phasedom),1:length(oridom),1:length(sfdom),1:length(colordom));

seqdum = [];
for l = 1:ceil(N_Im/length(immesh(:)))
    s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed+(l-1)*100);
    [dumval dumid] = sort(rand(s,1,length(immesh(:))));
    seqdum = [seqdum dumid];
end

seqdum = seqdum(1:N_Im);
imseq = immesh(seqdum);


%% Replace the last values of the sequence with blanks and then reshuffle.
%The last gratings of the sequence are the ones most likely to have multiple
%presentations.
blankflag = zeros(1,N_Im);
nblanks = round(P.blankProb*N_Im);

imseq((end-nblanks+1):end) = length(imdom) + 1;

s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed);
[dumval dumid] = sort(rand(s,1,length(imseq)));
imseq = imseq(dumid);    

blankflag(find(imseq == length(imdom)+1)) = 1;


%% Make the mask

%First get a sample image to id the dimensions
imgPath = fullfile(stimulusFolder, imdom(1).name);
Im = imread(imgPath);

[yN xN] = size(Im);
truncyN = round(size(Im,1)*P.truncYpercent/100);
truncxN = round(size(Im,2)*P.truncXpercent/100);
yN = truncyN*P.NyTile;
xN = truncxN*P.NxTile;
xdom = linspace(-P.x_size/2,P.x_size/2,xN);
ydom = linspace(-P.y_size/2,P.y_size/2,yN);
[xdom ydom] = meshgrid(xdom,ydom);
r = sqrt(xdom.^2 + ydom.^2);
if strcmp(P.mask_type,'disc')
    mask = zeros(size(r));
    id = find(r<=P.mask_radius);
    mask(id) = 1;
elseif strcmp(P.mask_type,'gauss')
    mask = exp((-r.^2)/(2*P.mask_radius^2));
elseif strcmp(P.mask_type,'cos')
    mask = cos(r/P.mask_radius*pi)/2 + .5;
    id = find(r>=P.mask_radius);
    mask(id) = 0;
else
    mask = [];
end
mask = single(mask);

%% Make smoothing kernel


% 
% if P.LPfilt    
%     
%     P.LPcutoff = 8; % half max: cyc/deg
%     LPradius = 1/(2*P.LPcutoff); %half cycle of raised cosine filter
%     
%     %First get a sample image to id the dimensions
%     imgPath = fullfile(stimulusFolder, imdom(1).name);
%     Im = imread(imgPath);
%     [yN xN] = size(Im);
%     
%     xdom = linspace(-P.x_size/2,P.x_size/2,xN);
%     ydom = linspace(-P.y_size/2,P.y_size/2,yN);
%     [xdom ydom] = meshgrid(xdom,ydom);
%     r = sqrt(xdom.^2 + ydom.^2);
%     
%     LPkern = cos(r/LPradius*pi)/2 + .5;
%     id = find(r>=LPradius);
%     LPkern(id) = 0;
%     
%    
%     LPkern = single(LPkern);
%     LPkern = LPkern/sum(LPkern(:));
%     LPkern = abs(fft2(LPkern));
%     
% end

nB = 4; %butterworth order

%or do it in the frequency domain
if P.LPfilt    
    
    %First get a sample image to id the dimensions
    %imgPath = fullfile(stimulusFolder, imdom(1).name);
    %Im = imread(imgPath);
    %[yN xN] = size(Im);
    xN = truncxN;
    yN = truncyN;
    xdom = linspace(0,xN/P.x_size,xN+1); xdom = xdom(1:end-1);
    xdom = xdom-xdom(length(xdom)/2+1);
    ydom = linspace(0,yN/P.y_size,yN+1); ydom = ydom(1:end-1);
    ydom = ydom-ydom(length(ydom)/2+1);
    [xdom ydom] = meshgrid(xdom,ydom);
    r = sqrt(xdom.^2 + ydom.^2);
    LPkern = 1./(1+(r/P.LPcutoff).^(2*nB));

    LPkern = fftshift(fftshift(LPkern,2),1); %put 0Hz in top left
    
%     LPkern = ifft2(LPkern);   
%     LPkern = LPkern/sum((LPkern(:)));
%     LPkern = fft2(LPkern);
end



%% Upload images to the card

nbit = 8; %hack. assumes loaded images are 8-bit.

for i = 1:N_Im    %loop through each image in the sequence
    if ~blankflag(i)
                   
        imgPath = fullfile(stimulusFolder, imdom(imseq(i)).name);
        Im = imread(imgPath);
        Im = single(Im); 
        
        Im = Im(1:truncyN,1:truncxN);
        
        if P.LPfilt
           Im = ifft2(fft2(Im).*LPkern); 
        end
        
        Im = 2*(Im/(2^nbit - 1) - .5); %-1 to 1... sets it up for ImtoRGB
        
        Im = Im*P.contrast/100;
        
        Im = repmat(Im,P.NyTile,P.NxTile);
        
        if P.distortbit
           Im = distortImage(Im,P); 
        end
        
        putinTexture(Im,P,i,mask); %Put in texture as RGB
    else
        putinTexture(0,P,i,mask); %Blank
    end
end
TDim = size(Im);

%Save it if 'running' experiment
if Mstate.running
    Pseq = struct;
    Pseq.imseq = imseq;
    
    domains = struct;
    domains.imdom = imdom;
    
    saveLog_savedImages(domains,Pseq)
    
    %This version of 'saveLog' comes from the old code, where I would call
    %it from 2 places (the make/play files).  It no longer needs the
    %varargin input structure, but I just kept it out of consistency (and laziness).
%     if loopTrial == 1
%         saveLog(domains)
%     else
%         saveLog(Pseq,P.rseed)  %append log file with the latest sequence
%     end
    
end



function putinTexture(Im,P,i,mask)

global Gtxtr screenPTR


Imcell{1} = Im;

Idraw = ImtoRGB(Imcell,P.colormod,P,mask);
Gtxtr(i) = Screen(screenPTR, 'MakeTexture', Idraw);
