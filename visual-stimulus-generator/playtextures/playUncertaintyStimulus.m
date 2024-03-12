function playUncertaintyStimulus

% yb, 06/15/2018, version 1.0

%This one uses the sequences that were already defined in the make file

global Mstate screenPTR screenNum daq loopTrial TimingInfo

global Gtxtr TDim  %Created in makeGratingTexture

global Stxtr %Created in makeSyncTexture

% yb: i feel bad doing this...(don't do this from now on)
global stimulusFolder scrambledFileNames scrambledContrasts stimOnset stimOffset 


P = getParamStruct;

screenRes       = Screen('Resolution',screenNum);
screenRes.hz    = Mstate.refresh_rate; %Replace it with this. The other sometimes reports 0hz.
pixpercmX       = screenRes.width/Mstate.screenXcm;
pixpercmY       = screenRes.height/Mstate.screenYcm;

syncWX          = round(pixpercmX*Mstate.syncSize);
syncWY          = round(pixpercmY*Mstate.syncSize);

white           = WhiteIndex(screenPTR); % pixel value for white
black           = BlackIndex(screenPTR); % pixel value for black
gray            = (white+black)/2;

%The following gives inaccurate spatial frequencies
% xN = 2*Mstate.screenDist*tan(P.x_size/2*pi/180);  %grating width in cm
% xN = round(xN*pixpercmX);  %grating width in pixels
% yN = 2*Mstate.screenDist*tan(P.y_size/2*pi/180);  %grating height in cm
% yN = round(yN*pixpercmY);  %grating height in pixels

%The following assumes the screen is curved
xcm     = 2*pi*Mstate.screenDist*P.x_size/360;  %stimulus width in cm
xN      = round(xcm*pixpercmX);  %stimulus width in pixels
ycm     = 2*pi*Mstate.screenDist*P.y_size/360;   %stimulus height in cm
yN      = round(ycm*pixpercmY);  %stimulus height in pixels

%Note: I used to truncate these things to the screen size, but it is not
%needed.  It also messes things up.
xran    = [P.x_pos-floor(xN/2)+1  P.x_pos+ceil(xN/2)];
yran    = [P.y_pos-floor(yN/2)+1  P.y_pos+ceil(yN/2)];

Npreframes  = ceil(P.predelay*screenRes.hz);
Npostframes = ceil(P.postdelay*screenRes.hz);


% N_Im = P.Nframes;
% 
% %Preallocate timing vectors
% n_frames_on = round(P.stim_on/(1000/Mstate.refresh_rate)); % convert ms to frames
% n_frames_off = round(P.stim_off/(1000/Mstate.refresh_rate)); % convert ms to frames
% 
% %Nflips = Npreframes + Npostframes + N_Im*P.h_per + N_Im*P.h_per_blanks + 1;
% Nflips = Npreframes + Npostframes + N_Im*n_frames_on + N_Im*n_frames_off + 1;
% VBLTimestamp = zeros(1,Nflips);
% StimulusOnsetTime = zeros(1,Nflips);
% FlipTimestamp = zeros(1,Nflips);

if loopTrial == 1
    TimingInfo = struct;
end

%%%%
%SyncLoc = [0 screenRes.height-syncWY syncWX-1 screenRes.height-1]';
SyncLoc = [0 0 syncWX-1 syncWY-1]';
SyncPiece = [0 0 syncWX-1 syncWY-1]';
StimLoc = [xran(1) yran(1) xran(2) yran(2)]';
srcrect = [0 0 TDim(1) TDim(2)]';


%Wake up the daq:
DaqDOut(daq, 0, 0); %I do this at the beginning because it improves timing on the first call to daq below

%f=1;

% save log: beforeOrAfter (0: before, 1: after)
beforeOrAfter = 0;
saveLog_uncertaintyStimulus(beforeOrAfter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%YB:  MAKE ON THE FLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nFiles                  = numel(scrambledFileNames);
BG_LUMINANCE            = 128;
normBG                  = BG_LUMINANCE/256;
        
%%% APERTURE MASK
% flat top 8
dummy                   = load(fullfile(stimulusFolder, scrambledFileNames{1}),'trial');
maskImgPath             = fullfile(stimulusFolder,'mask','Flattop8.tif');
maskImg                 = im2double(imread(maskImgPath)); % [0, 1]
flattop8                = imresize(squeeze(maskImg(:,:,end)), [512, 512]); % [0, 1]

% % save log
% saveLog_uncertaintyStimulus([]);
photodiode_on           = 1;
photodiode_off          = 2;

testTempGtxtr           = nan(20,1);
    

for iTrial = 1:nFiles
    % contrast
    iContrast = scrambledContrasts(iTrial); 
    
    % load variable 'trial'
    load(fullfile(stimulusFolder,scrambledFileNames{iTrial}), 'trial');
    
    
    [im_y, im_x, nFrames]   = size(trial);
    tempGtxtr               = nan(nFrames,1);
    
    
    trial = trial - normBG; % [-.5 .5]
    trial = trial .* iContrast;
    trial = trial .* repmat(flattop8, [1,1,nFrames]);
    
    for iFrame = 1:nFrames
        frameImg    = squeeze(trial(:,:,iFrame)); % [-0.5, 0.5]
        
        ImRGB       = zeros(im_y, im_x, 3,'uint8');  %make ImRGB uint8
        
        ImRdum      = frameImg;   %[-0.5 0.5]
        ImGdum      = frameImg;
        ImBdum      = frameImg;
       
        ImRdum = (ImRdum+normBG);  %[0 1]
        ImGdum = (ImGdum+normBG);
        ImBdum = (ImBdum+normBG);
        
        ImRGB(:,:,1)        = round(ImRdum*255);  %[0 255]
        ImRGB(:,:,2)        = round(ImGdum*255);
        ImRGB(:,:,3)        = round(ImBdum*255);
        tempGtxtr(iFrame)   = Screen(screenPTR, 'MakeTexture', ImRGB);
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% debugging (9/10/18)
    if(iTrial == 1)
        testTempGtxtr   = tempGtxtr;
    else
        tempGtxtr       = testTempGtxtr;
    end
    
    
    Screen(screenPTR, 'FillRect', P.background)
    
    %%%Play predelay %%%%
    %Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
    Screen('DrawTexture', screenPTR, Stxtr(photodiode_off),SyncPiece,SyncLoc);
    %[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    [VBLTimestamp, StimulusOnsetTime, FlipTimestamp] = Screen(screenPTR, 'Flip');
    %f=f+1;
    startTime = FlipTimestamp;
    if loopTrial ~= -1
        digWord = 7;  %Make 1st,2nd,3rd bits high
        DaqDOut(daq, 0, digWord);
    end
    for i = 2:Npreframes
        %Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
        Screen('DrawTexture', screenPTR, Stxtr(photodiode_off),SyncPiece,SyncLoc);
        Screen(screenPTR, 'Flip');
        %[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
        %f=f+1;
    end
    
    %%%%%Play whats in the buffer (the stimulus)%%%%%%%%%%
    
    %Unlike periodic grater, this doesn't produce a digital sync on last frame, just
    %the start of each grating.  But this one will always show 'h_per' frames on
    %the last grating, regardless of 'stimtime'.
    
    frame_index     = 1;
%     Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%     [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%     stimOnset(iTrial) = VBLTimestamp(f);
%     f=f+1;
%     Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%     [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%     f=f+1;
%     Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%     [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%     f=f+1;
    Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
    [VBLTimestamp, StimulusOnsetTime, FlipTimestamp] = Screen(screenPTR, 'Flip');
    stimOnset(iTrial) = FlipTimestamp;
    Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
    Screen(screenPTR, 'Flip');
    Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
    Screen(screenPTR, 'Flip');

        %digWord = bitxor(digWord,4);  %toggle only the 3rd bit on each grating update
        %DaqDOut(daq,0,digWord);
        
    for frame_index = 2:nFrames-1                  %sync flips on each update

%         Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;
%         
%         Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;
%         
%         Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;

        Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
        Screen(screenPTR, 'Flip');
        Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
        Screen(screenPTR, 'Flip');
        Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
        Screen(screenPTR, 'Flip');
    end
        
%         frame_index     = nFrames;
%         Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;
%         Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;
%         Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_off)],[],[StimLoc SyncLoc]);
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;

    frame_index     = nFrames;
    Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
    Screen(screenPTR, 'Flip');
    Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_on)],[],[StimLoc SyncLoc]);
    Screen(screenPTR, 'Flip');
    Screen('DrawTextures', screenPTR, [tempGtxtr(frame_index) Stxtr(photodiode_off)],[],[StimLoc SyncLoc]);
    [VBLTimestamp, StimulusOnsetTime, FlipTimestamp] = Screen(screenPTR, 'Flip');
    stimOffset(iTrial) = FlipTimestamp;
        %digWord = bitxor(digWord,4);  %toggle only the 3rd bit on each grating update
        %DaqDOut(daq,0,digWord);
        
        
%         %for j = 2:P.h_per_blanks                  %sync flips on each update
%         for j = 2:n_frames_off                  %sync flips on each update
%             Screen('DrawTextures', screenPTR, [Gtxtr(end) Stxtr(2-rem(i,2))],...
%                 [],[StimLoc SyncLoc]);
%             
%             [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%             f=f+1;
%         end
        
        
        
    %end
    
    
    
    %%%Play postdelay %%%%
    for i = 1:Npostframes-1
        %Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
        Screen('DrawTexture', screenPTR, Stxtr(photodiode_off),SyncPiece,SyncLoc);
        Screen(screenPTR, 'Flip');
%         [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
%         f=f+1;
    end
%     Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
    Screen('DrawTexture', screenPTR, Stxtr(photodiode_off),SyncPiece,SyncLoc);
    Screen(screenPTR, 'Flip');
    %[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    %f=f+1;
    
    %digWord = bitxor(digWord,7); %toggle all 3 bits (1st/2nd bits go low, 3rd bit is flipped)
    %DaqDOut(daq, 0,digWord);
    
    if loopTrial ~= -1
        DaqDOut(daq, 0, 0);  %Make sure 3rd bit finishes low
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    Screen('DrawTexture', screenPTR, Stxtr(photodiode_off),SyncPiece,SyncLoc);
    Screen(screenPTR, 'Flip');
    %[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    
%     %%Store timing information
%     if loopTrial ~= -1
%         TimingInfo.VBLTimestamp{loopTrial} = VBLTimestamp;
%         TimingInfo.StimulusOnsetTime{loopTrial} = StimulusOnsetTime;
%         TimingInfo.FlipTimestamp{loopTrial} = FlipTimestamp;
%         
%         %saveLog_uncertaintyStimulus([])
%     end
    
end

% save log: beforeOrAfter (0: before, 1: after)
beforeOrAfter = 1;
saveLog_uncertaintyStimulus(beforeOrAfter);

