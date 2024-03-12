function Displaycb(obj,event)
%Callback function from Stimulus PC.  Happens after 'Go' stimulus is
%executed.  Not after 'Build' stimulus.

global DcomState Mstate GUIhandles

n=get(DcomState.serialPortHandle,'BytesAvailable');
if n > 0
    inString = fread(DcomState.serialPortHandle,n);
    inString = char(inString');
else
    return
end

inString = inString(1:end-1)  %Get rid of the terminator

%'nextT' is the string sent after stimulus is played

%run2 should only be called here if it is not called within run2. It is
%called in run 2 if it needs to wait for something to happen between
%trials, such as saving or analysis.

% if strcmp(inString,'nextT') && Mstate.Ephys & get(GUIhandles.main.F0analysisFlag,'value')  %Why don't I do this in run2?  I think it had something to do with trial incrementing
%     electrode = Mstate.electrode;  %'A2x16'; %'A4x8'
%     %This will enter a while loop to wait
%     %for second TTL (trial end) before analysis
%     ephysOnlinePlotter(electrode)
% end

%run2 is not called here if in WF because it needs to wait for it to save
%the data between trials

A = Mstate.WF;
B = Mstate.Ephys;
C = get(GUIhandles.main.F0analysisFlag,'value');



if Mstate.running
    
    if strcmp(inString,'nextT')
        
        if B && C
            electrode = Mstate.electrode;  %'A2x16'; %'A4x8
            ephysOnlinePlotter(electrode) %Was having problems with the while loop that waits for TTL, so did this here instead of in run2
        end
        
        if ~A %Conditions when we don't need to wait between trials and should just let the signal from stimulator initiate the next trial.
            run2  %only run4 is called during widefield, which disables this callback function. run2 is called for everything else
            
            %run4
        end
    end
    
    % if strcmp(inString,'nextT') && Mstate.Ephys & get(GUIhandles.main.F0analysisFlag,'value')  %Why don't I do this in run2?  I think it had something to do with trial incrementing
    %     electrode = Mstate.electrode;  %'A2x16'; %'A4x8'
    %     %This will enter a while loop to wait
    %     %for second TTL (trial end) before analysis
    %     ephysOnlinePlotter(electrode)
    % end
end
