function makeUncertaintyStimulusTexture

% yb, 06/15/2018, version 1.0

%This one builds only the images that get played for this trial.
%Importantly, this one doesn't allow for drift or contrast reversal like
%the old version in the 'archive' file.

global Mstate screenPTR screenNum loopTrial % stimulus Block

global Gtxtr TDim%'playgrating' will use these

% yb: i feel bad doing this...(don't do this from now on)
global stimulusFolder scrambledFileNames scrambledContrasts stimOnset stimOffset

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


%% READ FROM STIMULUS FOLDER sca
% yb, 06/15/2018
stimulusFolder  = P.image_folder; %'/Users/iannauhaus/Desktop/uncertainty_stimulus';

% iterate through stimulus directory and save all filenames 
fileStructs             = dir(fullfile(stimulusFolder, 'D*.mat'));
fileNames               = {fileStructs.name};

% testing...
DEBUG_FLAG              = false;
if(DEBUG_FLAG)
    fileNames           = fileNames(1:3); % 16 seconds
    %fileNames           = fileNames(1:200); % 20 min
end

numFiles                = numel(fileNames);

scrambledFileNames      = [];
scrambledContrasts      = [];

% modulate contrasts (1, 0.38, 0.13)
contrastLevels          = [1, 0.39, 0.15, 0];
% testing...
if(DEBUG_FLAG)
    numBlanks           = 1;
else
    numBlanks           = 400; 
end

for i = 1:length(contrastLevels)
    if(contrastLevels(i) > 0)
        scrambledFileNames  = cat(2, scrambledFileNames, fileNames);
        scrambledContrasts  = cat(2, scrambledContrasts, repmat(contrastLevels(i),[1,numFiles]));
    else
        
        blankIndex          = 1:numBlanks;
        scrambledFileNames  = cat(2, scrambledFileNames, fileNames(blankIndex));
        scrambledContrasts  = cat(2, scrambledContrasts, repmat(contrastLevels(i),[1,numBlanks]));
    end
end

numTotalTrials          = numel(scrambledFileNames);
scrambleFileNameIndex   = randsample(numTotalTrials, numTotalTrials, false);
scrambleContrastIndex   = randsample(numTotalTrials, numTotalTrials, false);

scrambledFileNames      = scrambledFileNames(scrambleFileNameIndex);
scrambledContrasts      = scrambledContrasts(scrambleContrastIndex);  
stimOnset               = nan(numTotalTrials, 1);
stimOffset              = nan(numTotalTrials, 1);

dummy                   = load(fullfile(stimulusFolder, fileNames{1}),'trial');
[im_y, im_x, nFrames]   = size(dummy);

TDim = [im_y, im_x, nFrames];

% if Mstate.running
%     
%     Pseq = struct('imageID',fileNames,...
%                   'presentationIndex',shuffledIndex, ...
%                   'rseed',s, ...
%                   'userInputs',P);
%     
%     saveLog_naturalMovie(Pseq)
%     
%     %This version of 'saveLog' comes from the old code, where I would call
%     %it from 2 places (the make/play files).  It no longer needs the
%     %varargin input structure, but I just kept it out of consistency (and laziness).
% %     if loopTrial == 1
% %         saveLog(domains)
% %     else
% %         saveLog(Pseq,P.rseed)  %append log file with the latest sequence
% %     end
%     
% end


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
