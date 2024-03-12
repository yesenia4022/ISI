function configDisplayCom

global DcomState Mstate

%Modification of MP285Config, for configuration of udp port connection to visual stimulus PC  	

%DcomState is not initialized in the .ini file, nor is it saved in the state.headerString 

%close all open udp port objects on the same port and remove
%the relevant object form the workspace
port = instrfindall('RemoteHost',Mstate.stimulusUDP);
%port = instrfindall('RemoteHost','10.1.50.224');
% port = instrfindall('RemoteHost','128.83.22.157');

%port = '128.83.22.157';
if length(port) > 0; 
    fclose(port); 
    delete(port);
    clear port;
end

% make udp object named 'stim'
DcomState.serialPortHandle = udp(Mstate.stimulusUDP,'RemotePort',8866,'LocalPort',8844);
%DcomState.serialPortHandle = udp('10.1.50.224','RemotePort',8866,'LocalPort',8844);
% DcomState.serialPortHandle = udp('128.83.22.157','RemotePort',8866,'LocalPort',8844);


set(DcomState.serialPortHandle, 'OutputBufferSize', 1024)
set(DcomState.serialPortHandle, 'InputBufferSize', 1024)
set(DcomState.serialPortHandle, 'Datagramterminatemode', 'off')

%Establish serial port event callback criterion  
DcomState.serialPortHandle.BytesAvailableFcnMode = 'Terminator';
DcomState.serialPortHandle.Terminator = '~'; %Magic number to identify request from Stimulus ('c' as a string)
DcomState.serialPortHandle.bytesavailablefcn = @Displaycb;  

% open and check status 
fopen(DcomState.serialPortHandle);
stat=get(DcomState.serialPortHandle, 'Status');
if ~strcmp(stat, 'open')
    disp([' StimConfig: trouble opening port; cannot proceed']);
    DcomState.serialPortHandle=[];
    out=1;
    return;
end
