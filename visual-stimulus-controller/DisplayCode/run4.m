function run4

global GUIhandles Pstate Mstate looperInfo trialno syncInfo ACQ Tens tcourseSum TwoPcomState RipStruct

global DcomState

gigEflag = 0;

if Mstate.running %otherwise 'getnotrials' won't be defined for play sample
    nt = getnotrials;
end

if isfield(ACQ,'offLineSpatialBinning')
    D = ACQ.offLineSpatialBinning;
else
    D = 1;
end

if Mstate.running && trialno<=nt  %'trialno<nt' may be redundant.
    
    set(GUIhandles.main.showTrial,'string',['Trial ' num2str(trialno) ' of ' num2str(nt)] ), drawnow
    
    [c r] = getcondrep(trialno);  %get cond and rep for this trialno
    
    %%%Organization of commands is important for timing in this part of loop
    
    %tic
    if trialno == 1
        tic
        disp('Building stimulus...')
        buildStimulus(c,trialno)    %Tell stimulus to buffer the images (and move eye shutters)
        %Wait for serial port to respond from display. Make sure its
        %clear before waiting for bits to arrive next (leave input
        %blank)...
        waitforDisplayResp
        disp('Done building stimulus')
        buildstimulusTime = toc
        
    else %This new implementation pulls/saves/analyzes data of preveious trial while stimulus is building to save time
        disp('Building stimulus...')
        buildStimulus(c,trialno)    %Tell stimulus to buffer the images (and move eye shutters)
        
        %This safety measure clears the buffer before waiting for bits
        %to arrive. We want to do this right after build stimulus (i.e.
        %before display has a chance to respond) so that any bits
        %arriving in the buffer (in 'waitforDisplayResp_noclear below) are from
        %the display on this trial.
        
        clearDisplayBuffer
        %Why does data come into the display line after this?
        
        n = get(DcomState.serialPortHandle,'BytesAvailable')
                
        tic
        grabtimes = pullWFData(D);
        
        saveWFData(grabtimes,trialno-1)
        
        analyzeWFData(grabtimes,trialno-1)
        
        pull_Save_Analyze_Time = toc
        
        %Wait for bits to arrive on serial port from display when
        %stimulus is done building. It may already be done by the time
        %we get here. don't clear buffer.  Use '1' as input.
        waitforDisplayResp(1)
    end
    
    %ITDelay = toc;
    %disp(['Intertrial Delay = ' num2str(ITDelay)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if ACQ.btwTrialShutter
        str = 'mm\r'; %open shutter
        LUMsend(str)
    end
    if ACQ.MechbtwTrialShutter
        writePosition(ACQ.servo,ACQ.servoOpenVal)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    startStimulus      %Tell Display to show its buffered images. TTL from stimulus computer "feeds back" to frame grabber
    disp('Playing stimulus')
        
    
    pause(.5) %Make sure it has had time to receive TTL from vstim, otherwise it won't enter "logging" state in time
    
    display('logging data')
    
    
    try
        tic
        wait(ACQ.vid,ACQ.total_time+10,'logging')
        toc
        display('done logging data')
        
        %Need to wait for Display bits!!  Displaycb will not get called
        %while we are running in this loop. There were leftover bits
        %screwing things up on the next trial.
        %waitforDisplayResp(1)
        'hello1'
        pause(.5)
        %This is a hack to re-establish communication
        fwrite(DcomState.serialPortHandle,'X;dum;~')
        waitforDisplayResp(1)
        'hello2'
    catch
        ACQ.vid.FramesAvailable
        'Did not finish logging. Did it get all the frames?'
    end
    
    if gigEflag
        stop(ACQ.vid)
    end
    
    
    if ACQ.btwTrialShutter
        str = 'zz\r'; %close shutter
        LUMsend(str)
    end
    if ACQ.MechbtwTrialShutter
        writePosition(ACQ.servo,ACQ.servoCloseVal)
    end
    
    %If its the last trial, do all this time-consuming stuff now.
    %Otherwise do it while the stimulus is getting built above.
    %1) Pull data in.
    if nt == trialno
        
        grabtimes = pullWFData(D);
        saveWFData(grabtimes,trialno)
        analyzeWFData(grabtimes,trialno)
        
    end
    
    trialno = trialno+1;
    
    %If nothing needs to be done between trials (saving or analysis), run2 gets called at the end of the stimulus by
    %Displaycb.m. Otherwise...
        
    run4
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %Happens at the end of the experiment if not running wide-field acq.
    
    %set(GUIhandles.param.playSample,'enable','off')
    
    Mstate.running = 0;
    set(GUIhandles.main.runbutton,'string','Run')
        
end


function y = LUMsend(str)

%Send stuff to the luminator
global LUMserial

if(~isempty(LUMserial))
    fwrite(LUMserial,[str 13]);
    pause(0.2);
    N = get(LUMserial,'BytesAvailable');
    if(N>0)
        y = fgetl(LUMserial);
    end
else
    y = '\n Error>> No message from luminator!\n';
end

function y = DSimage(x,D)
y = 0;
for ii = 1:D
    for jj = 1:D
        y = y+x(ii:D:end,jj:D:end,:)/(D^2);
    end
end



function grabtimes = pullWFData(D)

global ACQ Tens

disp('Pulling data into Matlab')
[imdum grabtimes] = getdata(ACQ.vid,ACQ.vid.FramesAvailable);
if isempty(imdum)
    error('Did not grab frames!')
end
Tens = DSimage(imdum,D);
disp('Done pulling data into Matlab')

function saveWFData(grabtimes,trialno)

global Tens Mstate syncInfo

disp(['Saving trial: ' num2str(size(Tens,3)) ' frames'])
SaveTrial(trialno) %Otherwise save right after starting the stimulus build at the top of this function
syncInfo.acqSyncs = grabtimes;
syncInfo.frameRate = Mstate.framerate;  %frame rate of slave
saveSyncInfo(syncInfo,trialno)  %append .analyzer file
disp('Done saving.')


function analyzeWFData(grabtimes,trialno)

global Tens

%Online analysis
disp('Online analysis')
onlineF1Analysis(Tens,grabtimes,trialno)     %Compute F1
onlineF0Analysis(Tens,trialno)     %Compute F0
plot_WF_Timecourse(trialno)
disp('Done with online analysis')