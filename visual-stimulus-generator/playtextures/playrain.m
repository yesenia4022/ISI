function playrain_beta


global Mstate screenPTR screenNum loopTrial TimingInfo

global GtxtrAll StimLoc Dim_source Dim_dest daq  %Created in makeGratingTexture

global Stxtr %Created in makeSyncTexture


P = getParamStruct;

screenRes = Screen('Resolution',screenNum);

screenRes.hz = Mstate.refresh_rate;
pixpercmX = screenRes.width/Mstate.screenXcm;
pixpercmY = screenRes.height/Mstate.screenYcm;

syncWX = round(pixpercmX*Mstate.syncSize);
syncWY = round(pixpercmY*Mstate.syncSize);
SyncLoc = [0 0 syncWX-1 syncWY-1]';
SyncPiece = [0 0 syncWX-1 syncWY-1]';

StimPiece = [0 0 Dim_source(2)-1 Dim_source(1)-1]';
StimPiece = StimPiece*ones(1,P.Ndrops);

screenRes = Screen('Resolution',screenNum);
Npreframes = ceil(P.predelay*screenRes.hz);
Npostframes = ceil(P.postdelay*screenRes.hz);

%Preallocate timing vectors
Nflips = Npreframes+Npostframes+length(GtxtrAll)*P.h_per+1;
VBLTimestamp = zeros(1,Nflips);
StimulusOnsetTime = zeros(1,Nflips);
FlipTimestamp = zeros(1,Nflips);

if loopTrial == 1
    TimingInfo = struct;
end

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

%Wake up the daq:
%DaqDOut(daq, 0, 0); %I do this at the beginning because it improves timing on the first call to daq below

f = 1;

Screen(screenPTR, 'FillRect', P.background)

%%%Play predelay %%%%
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
f = f+1;
if loopTrial ~= -1
    digWord = 7;  %Make 1st,2nd,3rd bits high
    %DaqDOut(daq, 0, digWord);
end
for i = 2:Npreframes
    Screen('DrawTexture', screenPTR, Stxtr(2), SyncPiece, SyncLoc);
    [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    f = f+1;
end

%%%%%Play whats in the buffer (the stimulus)%%%%%%%%%%

for i = 1:length(GtxtrAll)
    
    Screen('DrawTextures', screenPTR, [GtxtrAll{i} Stxtr(2-rem(i,2))],...
        [StimPiece SyncPiece],[StimLoc{i,1} SyncLoc]);  
    
    [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    f = f+1;
    
%     digWord = bitxor(digWord,4);  %toggle only the 3rd bit on each grating update
%     DaqDOut(daq,0,digWord); 
    
    for j = 2:P.h_per                  %sync flips on each update    
        
        Screen('DrawTextures', screenPTR, [GtxtrAll{i} Stxtr(2-rem(i,2))],...
            [StimPiece SyncPiece],[StimLoc{i,j} SyncLoc]);          
        [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
        f = f+1;
        
    end

end

%%%Play postdelay %%%%
for i = 1:Npostframes-1
    Screen('DrawTexture', screenPTR, Stxtr(2),SyncPiece,SyncLoc);
    [VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
    f = f+1;
end
Screen('DrawTexture', screenPTR, Stxtr(1),SyncPiece,SyncLoc);
[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');
f = f+1;

if loopTrial ~= -1
    digWord = bitxor(digWord,7); %toggle all 3 bits (1st/2nd bits go low, 3rd bit is flipped)
    %DaqDOut(daq, 0,digWord);
    %DaqDOut(daq, 0, 0);  %Make sure 3rd bit finishes low
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Screen('DrawTexture', screenPTR, Stxtr(2),[0 0 syncWX-1 syncWY-1],SyncLoc);  
[VBLTimestamp(f) StimulusOnsetTime(f) FlipTimestamp(f)] = Screen(screenPTR, 'Flip');

%%Store timing information
if loopTrial ~= -1
    TimingInfo.VBLTimestamp{loopTrial} = VBLTimestamp;
    TimingInfo.StimulusOnsetTime{loopTrial} = StimulusOnsetTime;
    TimingInfo.FlipTimestamp{loopTrial} = FlipTimestamp;
    
    saveLog_rain([],[])
end

