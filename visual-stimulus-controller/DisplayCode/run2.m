function run2

global GUIhandles Pstate Mstate looperInfo trialno syncInfo ACQ Tens tcourseSum TwoPcomState RipStruct

gigEflag = 0;
if isfield(ACQ,'camera')
    if strcmp(ACQ.camera,'AVT') || strcmp(ACQ.camera,'Panda')
        gigEflag = 1;
    end
end

if Mstate.running %otherwise 'getnotrials' won't be defined for play sample
    nt = getnotrials;
end

if Mstate.running && trialno<=nt  %'trialno<nt' may be redundant.
    
    set(GUIhandles.main.showTrial,'string',['Trial ' num2str(trialno) ' of ' num2str(nt)] ), drawnow
    
    [c r] = getcondrep(trialno);  %get cond and rep for this trialno
    
    %     if Mstate.WF
    %         start(analogIN)  %Start sampling acquistion and stimulus syncs
    %     end

    if Mstate.Ephys
        xippmex('digin'); %First clear the dig buffer
    end
    
    %%%Organization of commands is important for timing in this part of loop
    
    disp('Building stimulus...')
    tic
    buildStimulus(c,trialno)    %Tell stimulus to buffer the images (and move eye shutters)
    
    waitforDisplayResp   %Wait for serial port to respond from display
    
    disp('Done building stimulus')
    toc
    if Mstate.WF
        if ACQ.btwTrialShutter
            str = 'mm\r'; %open shutter
            LUMsend(str)
        end
        if ACQ.MechbtwTrialShutter
            writePosition(ACQ.servo,ACQ.servoOpenVal) 
        end
    end
  
    startStimulus      %Tell Display to show its buffered images. TTL from stimulus computer "feeds back" to frame grabber
    disp('Playing stimulus')
    
    if Mstate.WF
        
        pause(.5) %Make sure it has had time to receive TTL from vstim, otherwise it won't enter "logging" state in time
        tic
        display('logging data')
                 
        if gigEflag
            while islogging(ACQ.vid)

            end
        end
        
        try

            wait(ACQ.vid,ACQ.total_time+10,'logging')
            
            display('done logging data')
            
            %Should wait for Display bits.  Displaycb will not get called
            %while we are running in this loop. 
            %waitforDisplayResp(1)
            
        catch
            ACQ.vid.FramesAvailable        
            'Did not finish logging. Did it get all the frames?'
        end        

        if gigEflag
            stop(ACQ.vid)
        end
        
        if ACQ.vid.FramesAvailable<ACQ.vid.FramesperTrigger
            display('Did not collect all the requested frames')
        end
        
        display('done logging data')
        
        if ACQ.btwTrialShutter
            str = 'zz\r'; %close shutter
            LUMsend(str)
        end
        if ACQ.MechbtwTrialShutter
            writePosition(ACQ.servo,ACQ.servoCloseVal)
        end       
        
        if isfield(ACQ,'offLineSpatialBinning')
            D = ACQ.offLineSpatialBinning;
        else 
            D = 1;
        end
        toc
        
        tic
        disp('Pulling data into Matlab')
        [grab_raw grabtimes] = getdata(ACQ.vid,ACQ.vid.FramesAvailable);
        
        if isempty(grab_raw)
            error('Did not grab frames!')
        end
        
        Tens = DSimage(grab_raw,D);
        disp('Done pulling data into Matlab')
        
        toc
        
        if gigEflag
            flushdata(ACQ.vid) % Need to do this for gigE camera case
            updateTempReading
            
            start(ACQ.vid) %Arm the hardware trigger
        end
        
        
        if ACQ.timecourseBit
            if r == 1
                tcourseSum{c} = 0;
            end
            superPixID = find(ACQ.timecourseROI);
            for q = 1:size(Tens,3)
                imdum = Tens(:,:,q);
                tcourse(q) = mean(imdum(superPixID));
            end
            tcourseSum{c} = tcourseSum{c} + tcourse;
            
            figure(99),
            
            nc = getnoconds;
            tdom = (0:length(tcourse)-1)/ACQ.FPS;
            subplot(2,nc,c), 
            %subplot(ceil(sqrt(nc)),ceil(sqrt(nc)),c), 
            plot(tdom,tcourse(1:end)), hold on            
            plot(tdom,tcourseSum{c}/r,'r'), xlabel('sec'), hold off 
            if c ==1
            legend('trial','average')
            end
            title(['condition' num2str(c)])
            
            subplot(2,nc,c+nc), imagesc(squeeze(mean(squeeze(Tens(:,:,2:end-1)),3))), axis image
            

% figure(99), imagesc(squeeze(mean(squeeze(Tens(:,:,2:end-1)),3))), axis image, colormap gray
% figure(100), plot(tcourse)
        end


        %%
        disp(['Saving trial: ' num2str(size(Tens,3)) ' frames'])
tic
        SaveTrial(trialno)

        disp('Done saving.')
        
        syncInfo.acqSyncs = grabtimes;
        syncInfo.frameRate = Mstate.framerate;  %frame rate of slave
        saveSyncInfo(syncInfo,trialno)  %append .analyzer file
      
%        savetime = toc
toc
    end
    
    %Tens is a global variable to give it a dedicated memory allocation.
    
    if Mstate.WF     
        onlineF1Analysis(Tens,grabtimes,trialno)     %Compute F1        
        onlineF0Analysis(Tens,trialno)     %Compute F0        
    end
    
    
    trialno = trialno+1;
    
    %If nothing needs to be done between trials (saving or analysis), run2 gets called at the end of the stimulus by
    %Displaycb.m. Otherwise...

    if Mstate.WF
        run2 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %Happens at the end of the experiment if not running wide-field acq.
   
    
    %set(GUIhandles.param.playSample,'enable','off')
    
    Mstate.running = 0;
    set(GUIhandles.main.runbutton,'string','Run')
    
    if Mstate.twoP
        fprintf(TwoPcomState.serialPortHandle,'S');  %Stop acquiring
    end
    
    if Mstate.Ephys
        
        status = xippmex;
        RipStruct
        if status ~= 1;
            error('Xippmex Did Not Initialize');
        end
        %[trial_descriptor] = xippmex('trial', RipStruct.opers(1), 'stopped')
        [trial_descriptor] = xippmex('trial','stopped')
        %xippmex('close');
        
       
        %         %[trial_descriptor] = xippmex('trial', RipStruct.opers, 'stopped', RipStruct.Filename);
        %          RipStruct.opers
        %          opers = xippmex('opers');
        
        %          %xippmex('trial', opers(1), 'stopped');
        %          if ~isempty(opers)
        %              desc = xippmex('trial', opers(1));
        %              while ~strcmp(desc.status, 'stopped')
        %                  pause(1); desc = xippmex('trial', opers(1));
        %                  xippmex('trial', RipStruct.opers(1), 'stopped');
        %              end
        %              xippmex('trial', opers(1), 'stopped');
        %          else
        %              pause(1)
        %              xippmex('close');
        %          end
        
    end
    
    
%     if get(GUIhandles.main.widefieldflag,'value')
%         %set(GUIhandles.param.playSample,'enable','off')
%         saveOnlineAnalysis
%     end
    
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