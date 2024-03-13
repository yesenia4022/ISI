function playnaturalmovie

%This one uses the sequences that were already defined in the make file

global Mstate screenPTR screenNum daq loopTrial TimingInfo

global Gtxtr TDim  %Created in makeGratingTexture

global Stxtr %Created in makeSyncTexture

global imageStruct imageOnset

P = getParamStruct;

screenRes = Screen('Resolution',screenNum);
screenRes.hz = Mstate.refresh_rate; %Replace it with this. The other sometimes reports 0hz.
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

syncWX = round(pixpercmX*Mstate.syncSize);
syncWY = round(pixpercmY*Mstate.syncSize);


white = WhiteIndex(screenPTR); % pixel value for white
black = BlackIndex(screenPTR); % pixel value for black
gray = (white+black)/2;

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

%Note: I used to truncate these things to the screen size, but it is not
%needed.  It also messes things up.
xran = [P.x_pos-floor(xN/2)+1  P.x_pos+ceil(xN/2)];
yran = [P.y_pos-floor(yN/2)+1  P.y_pos+ceil(yN/2)];

Npreframes = ceil(P.predelay*screenRes.hz);
Npostframes = ceil(P.postdelay*screenRes.hz);

N_Im = P.Nframes;

%Preallocate timing vectors
n_frames_on = round(P.stim_on/(1000/Mstate.refresh_rate)); % convert ms to frames
n_frames_off = round(P.stim_off/(1000/Mstate.refresh_rate)); % convert ms to frames
    
%Nflips = Npreframes + Npostframes + N_Im*P.h_per + N_Im*P.h_per_blanks + 1;
Nflips = Npreframes + Npostframes + N_Im*n_frames_on + N_Im*n_frames_off + 1;
VBLTimestamp = zeros(1,Nflips);
StimulusOnsetTime = zeros(1,Nflips);
FlipTimestamp = zeros(1,Nflips);

if loopTrial == 1
    TimingInfo = struct;
end

%%%%
%SyncLoc = [0 screenRes.height-syncWY syncWX-1 screenRes.height-1]';
SyncLoc = [0 0 syncWX-1 syncWY-1]';
SyncPiece = [0 0 syncWX-1 syncWY-1]';
StimLoc = [xran(1) yran(1) xran(2) yran(2)]';
StimPiece = [0 0 TDim(1)-1 TDim(2)-1]';
% srcrect = [0 0 TDim(1) TDim(2)]';


%Wake up the daq:
%DaqDOut(daq, 0, 0); %I do this at the beginning because it improves timing on the first call to daq below

f=1;

Screen(screenPTR, 'FillRect', P.background)

%%%Play predelay %%%%
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
f=f+1;
if loopTrial ~= -1
    digWord = 7;  %Make 1st,2nd,3rd bits high
    %DaqDOut(daq, 0, digWord);
end
for i = 2:Npreframes
    Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    f=f+1;
end

%%%%%Play whats in the buffer (the stimulus)%%%%%%%%%%

%Unlike periodic grater, this doesn't produce a digital sync on last frame, just
%the start of each grating.  But this one will always show 'h_per' frames on
%the last grating, regardless of 'stimtime'.
    
%for i = 1:N_Im
for i = 1 : length(Gtxtr)
    
    for j = 1
        Screen('DrawTextures', screenPTR, [Gtxtr(i) Stxtr(1)],[StimPiece SyncPiece],[StimLoc SyncLoc]);
        [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
        imageOnset(i) = VBLTimestamp(f);
        f=f+1;
    end
    %digWord = bitxor(digWord,4);  %toggle only the 3rd bit on each grating update
    %DaqDOut(daq,0,digWord);
    
    
    %for j = 2:P.h_per                  %sync flips on each update
    for j = 2:n_frames_on                  %sync flips on each update
        Screen('DrawTextures', screenPTR, [Gtxtr(i) Stxtr(2)],...
            [StimPiece SyncPiece],[StimLoc SyncLoc]);

        [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
        f=f+1;
    end
    
    for j = 1
        Screen('DrawTextures', screenPTR, [Gtxtr(end) Stxtr(2)],[StimPiece SyncPiece],[StimLoc SyncLoc]);
        [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
        f=f+1;
    end
    %digWord = bitxor(digWord,4);  %toggle only the 3rd bit on each grating update
    %DaqDOut(daq,0,digWord);
    
    
    %for j = 2:P.h_per_blanks                  %sync flips on each update
    for j = 2:n_frames_off                  %sync flips on each update
        Screen('DrawTextures', screenPTR, [Gtxtr(end) Stxtr(2)],...
            [StimPiece SyncPiece],[StimLoc SyncLoc]);
        [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
        f=f+1;

    end
    
    
    
end

    

%%%Play postdelay %%%%
for i = 1:Npostframes-1
    Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    f=f+1;
end
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
f=f+1;
%digWord = bitxor(digWord,7); %toggle all 3 bits (1st/2nd bits go low, 3rd bit is flipped)
%DaqDOut(daq, 0,digWord);  

if loopTrial ~= -1
    %DaqDOut(daq, 0, 0);  %Make sure 3rd bit finishes low
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);  
[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');

%%Store timing information
if loopTrial ~= -1
    TimingInfo.VBLTimestamp{loopTrial} = VBLTimestamp;
    TimingInfo.StimulusOnsetTime{loopTrial} = StimulusOnsetTime;
    TimingInfo.FlipTimestamp{loopTrial} = FlipTimestamp;
    
    saveLog_naturalMovie(imageStruct)
end

