function makeSavedImageTexture

%This one builds only the images that get played for this trial.
%Importantly, this one doesn't allow for drift or contrast reversal like
%the old version in the 'archive' file.

global Mstate screenNum loopTrial %movieBlock

global Gtxtr TDim  %'playgrating' will use these


P = getParamStruct;

if loopTrial <= 1 || P.repeatTrial_bit == 0
    Screen('Close')  %First clean up: Get rid of all textures/offscreen windows   
    Gtxtr = []; TDim = [];  %reset
end

screenRes = Screen('Resolution',screenNum); %sometimes this reports 0Hz
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

%Replace it with this (comes from flipInterval)
screenRes.hz = Mstate.refresh_rate;



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Create Sequence for this trial%%%%%

%colordom = getColorDomain(P.colorspace);

%number of images to present in entire trial
N_Im = round(P.stim_time*screenRes.hz/P.h_per); 


%%
stimulusFolder = P.image_folder; %'/Users/iannauhaus/Desktop/Textures';

trialListA = dir(fullfile(stimulusFolder, P.fileID)); %Conditioned list in folder
trialListB = dir(fullfile(stimulusFolder, P.fileID_block2)); %Conditioned list in folder

%Implemented 11/29/20
dumListA = natsort({trialListA.name});
dumListB = natsort({trialListB.name});
for i = 1:length(dumListA)
    trialListA(i).name = dumListA{i};
    trialListB(i).name = dumListB{i}; 
end
%%%%%%%%%%%%%%%%%%%%%%

if ~(P.block_bit)
    imdom = trialListA;  %image "domain" is the list of potential image file names
else
    imdom = [trialListA; trialListB]; %If we want to run alternating blocks [a b a b...]
end

%%
immesh = ndgrid(1:length(imdom)); %degenerate use of ndgrid, but keeping consistent with other modules
%[phasemesh orimesh sfmesh colormesh] = ndgrid(1:length(phasedom),1:length(oridom),1:length(sfdom),1:length(colordom));


if ~(P.block_bit)
    block_Nimage = N_Im; %if we are not running blocks, they are the same.
else
    block_Nimage = P.block_Nimage;  
end

blockseq = [zeros(1,block_Nimage) ones(1,block_Nimage)];
blockseq = repmat(blockseq,[1 ceil(N_Im/length(blockseq))]);
blockseq = blockseq(1:N_Im);

for i = 1:N_Im %Loop each flashed image (presented h_per frames)
    
    %Used this until 8/10/18
    %s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed+i*20);
    
    if P.rseed > 0
        
        s = RandStream('mt19937ar','Seed',P.rseed*100 + i*20);
        
        if ~blockseq(i)
            idx = round(rand(s,1,1)*(length(trialListA)-1)+1); %take a single random sample and "replace"
        else
            idx = round(rand(s,1,1)*(length(trialListB)-1)+1);
            idx = idx+length(trialListA);
        end
        
    else  %don't randomize
        
        idx = rem(i,length(imdom));  %If there are more images played than are in the folder       
        if idx == 0  
            idx = length(imdom);           
        end
        
    end
    
    seqdum(i) = idx;
        
end

imseq = immesh(seqdum);


immesh = ndgrid(1:length(imdom)); %degenerate use of ndgrid, but keeping consistent with other modules
%[phasemesh orimesh sfmesh colormesh] = ndgrid(1:length(phasedom),1:length(oridom),1:length(sfdom),1:length(colordom));


%% randomly replace images with a blank

blankflag = zeros(1,N_Im); 
% nblanks = round(P.blankProb*N_Im);
% 
% imseq((end-nblanks+1):end) = length(imdom) + 1;
% 
% s = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed);
% [dumval dumid] = sort(rand(s,1,length(imseq)));
% imseq = imseq(dumid);    
% 
% blankflag(find(imseq == length(imdom)+1)) = 1;


%% Make the mask

%First get a sample image to id the dimensions
imgPath = fullfile(stimulusFolder, imdom(1).name);
Im = imread(imgPath);
Im = Im(:,:,1);

[yN xN] = size(Im);
truncyN = round(size(Im,1)*P.truncYpercent/100);
truncxN = round(size(Im,2)*P.truncXpercent/100);
%yN = truncyN*P.NyTile;
%xN = truncxN*P.NxTile;
yN = truncyN;
xN = truncxN;
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


if loopTrial <= 1 || P.repeatTrial_bit == 0
    
    nbit = 8; %hack. assumes loaded images are 8-bit.
    
    clear Im
    for i = 1:N_Im    %loop through each image in the sequence
        
        if ~blankflag(i)
            
            imgPath = fullfile(stimulusFolder, imdom(imseq(i)).name);
            Imdum = imread(imgPath);
            if P.gray_bit
                Imdum = mean(Imdum,3);
            end
            Imdum = single(Imdum);
            if length(size(Imdum)) == 3
                Im{1} = Imdum(:,:,1);
                Im{2} = Imdum(:,:,2);
                Im{3} = Imdum(:,:,3);
            else
                Im{1} = Imdum;
            end
            
            for cdim = 1:length(Im)
                
                %Im{cdim} = Im{cdim}(1:truncyN,1:truncxN);
                
                if P.LPfilt
                    Im{cdim} = ifft2(fft2(Im{cdim}).*LPkern);
                end
                
                Im{cdim} = 2*(Im{cdim}/(2^nbit - 1) - .5); %-1 to 1... sets it up for ImtoRGB
                
                Im{cdim} = Im{cdim}*P.contrast/100;
                
                %Im = repmat(Im,P.NyTile,P.NxTile);
                
                if P.distortbit
                    Im{cdim} = distortImage(Im{cdim},P);
                end
            end
            
            putinTexture(Im,P,i,mask); %Put in texture as RGB
        else
            putinTexture(0,P,i,mask); %Blank
        end
    end
    
    TDim = size(Im{1});
    
end

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



function putinTexture(Imcell,P,i,mask)

global Gtxtr screenPTR


Idraw = ImtoRGB(Imcell,P.colormod,P,mask);
Gtxtr(i) = Screen(screenPTR, 'MakeTexture', Idraw);
