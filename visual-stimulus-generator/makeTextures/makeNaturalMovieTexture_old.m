function makeNaturalMovieTexture

%This one builds only the images that get played for this trial.
%Importantly, this one doesn't allow for drift or contrast reversal like
%the old version in the 'archive' file.

global Mstate screenPTR screenNum loopTrial %movieBlock

global Gtxtr TDim%'playgrating' will use these

global stimOnset


Screen('Close')  %First clean up: Get rid of all textures/offscreen windows

Gtxtr = []; TDim = [];  %reset

P = getParamStruct;

screenRes = Screen('Resolution',screenNum); %sometimes this reports 0Hz
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

%Replace it with this (comes from flipInterval)
screenRes.hz = Mstate.refresh_rate;

%The following gives inaccurate spatial frequencies
% xN = 2*Mstate.screenDist*tan(P.x_size/2*pi/180);  %grating width in cm
% xN = round(xN*pixpercmX);  %grating width in pixels
% yN = 2*Mstate.screenDist*tan(P.y_size/2*pi/180);  %grating height in cm
% yN = round(yN*pixpercmY);  %grating height in pixels

%The following assumes the screen is curved
xcm = 2*pi*Mstate.screenDist*P.x_size/360;  %stimulus width in cm
xN = round(xcm*pixpercmX);  %stimulus width in pixels
ycm = 2*pi*Mstate.screenDist*P.y_size/360;   %stimulus height in cm
yN = round(ycm*pixpercmY);  %stimulus height in pixels


%% READ FROM STIMULUS FOLDER 
stimulusFolder  = P.image_folder; %'/Users/iannauhaus/Desktop/natural_stimulus';

% iterate through all images and save them into a matrix
imgIDs          = [];
Im              = [];
images          = dir(fullfile(stimulusFolder, '*.png'));
for i = 1:length(images)
%for i = 1:10
    % unique id
    imgIDs      = cat(1, imgIDs, str2double(images(i).name(1:3)));
    
    % file path
    imgPath   = fullfile(stimulusFolder, images(i).name);

    % save image to massive matrix
    Im          = cat(3, Im, imread(imgPath));
end


% shuffle images
N_Im            = size(imgIDs,1);
s               = RandStream.create('mrg32k3a','NumStreams',1,'Seed',P.rseed + loopTrial);
rng(s.Seed);
shuffledIndex   = Shuffle(1:N_Im)';

%% APERTURE MASK
% flat top 8
[im_y, im_x, ~] = size(Im);
maskImgFolder   = stimulusFolder;
maskImgPath     = fullfile(maskImgFolder,'mask','Flattop8.tif');
maskImg         = im2double(imread(maskImgPath));
maskImg         = imresize(squeeze(maskImg(:,:,end)), [im_y, im_x]);

% % raised-cosine
% [im_y, im_x, ~] = size(Im);
% [yy,xx]     = meshgrid((1:im_y)-im_y/2,(1:im_x)-im_x/2);
% radius      = ceil(min([im_y,im_x])/2);
% distMatrix  = sqrt(yy.^2 + xx.^2);
% distMatrix(distMatrix>radius) = radius;
% raisedCosine=(cos(distMatrix/radius * pi)+1)/2;

bitDepth                = 8;
pixRange                = 2^bitDepth - 1;
meanLum                 = 2^(bitDepth - 1);              
maxRangeZeroCentered    = 2^(bitDepth - 1); 
for i = 1:length(shuffledIndex)    %loop through each image in the sequence
    
    imgIndex        = shuffledIndex(i);
    
    normalizedIm    = double(squeeze(Im(:,:,imgIndex))) - P.background; % [0    255]
    normalizedIm    = normalizedIm ./ maxRangeZeroCentered;      % [-.5  +.5]
    normalizedIm    = normalizedIm .* maskImg;

    putinTexture(normalizedIm,P,i); %Put in texture as RGB
end
putinTexture(0,P,i+1); %Last element is always the Blank

stimOnset = nan(length(Gtxtr),1);

TDim = size(Im);

if Mstate.running
    
    Pseq = struct('imageID',imgIDs,...
                  'presentationIndex',shuffledIndex, ...
                  'rseed',s, ...
                  'userInputs',P);
    
    saveLog_naturalMovie(Pseq)
    
    %This version of 'saveLog' comes from the old code, where I would call
    %it from 2 places (the make/play files).  It no longer needs the
    %varargin input structure, but I just kept it out of consistency (and laziness).
%     if loopTrial == 1
%         saveLog(domains)
%     else
%         saveLog(Pseq,P.rseed)  %append log file with the latest sequence
%     end
    
end


function putinTexture(Im,P,i)

global Gtxtr screenPTR

Imcell{1} = Im;

% YB: Not using 'ImtoRGB'. Need to dynamically adjust overall luminance
% level (instead of fixing it to 0.5, or 128).
%Idraw = ImtoRGB(Imcell,1,P,[]);
ImRGB = zeros(length(Imcell{1}(:,1)),length(Imcell{1}(1,:)),3,'uint8');  %make ImRGB uint8

ImRdum = Imcell{1}*P.redgain;   %[-1 1]
ImGdum = Imcell{1}*P.greengain;
ImBdum = Imcell{1}*P.bluegain;

% if ~isempty(mask) %Made 'if' to significantly reduce computation for large stimuli
%     %C = (128 - P.background)/128;  %used so the mask looks right when background ~= 128
%     C = (P.background)/128;  %used so the mask looks right when background ~= P.background
%     ImRdum = (ImRdum + C).*mask - C;
%     ImGdum = (ImGdum + C).*mask - C;
%     ImBdum = (ImBdum + C).*mask - C;
% end
normBG = (P.background)/128;
% ImRdum = (ImRdum+1)/2 - (.5-P.redbase);  %[0 1]
% ImGdum = (ImGdum+1)/2 - (.5-P.greenbase);  
% ImBdum = (ImBdum+1)/2 - (.5-P.bluebase);  

ImRdum = (ImRdum+normBG)/2 - (.5-P.redbase);  %[0 1]
ImGdum = (ImGdum+normBG)/2 - (.5-P.greenbase);  
ImBdum = (ImBdum+normBG)/2 - (.5-P.bluebase);  

ImRGB(:,:,1) = round(ImRdum*255);  %[0 255]
ImRGB(:,:,2) = round(ImGdum*255);
ImRGB(:,:,3) = round(ImBdum*255);


Gtxtr(i) = Screen(screenPTR, 'MakeTexture', ImRGB);
