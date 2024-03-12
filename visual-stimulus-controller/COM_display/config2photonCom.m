function config2photonCom

global TwoPcomState Mstate

%Modification of MP285Config, for configuration of udp port connection to two photon PC (scanbox)  	

%TwoPcomState is not initialized in the .ini file, nor is it saved in the state.headerString 

%close all open udp port objects on the same port and remove
%the relevant object form the workspace
port = instrfindall('RemoteHost',Mstate.twoPhotonUDP);
%port = instrfindall('RemoteHost','10.1.50.224');
%port = '10.1.50.224';
if length(port) > 0; 
    fclose(port); 
    delete(port);
    clear port;
end

% make udp object
TwoPcomState.serialPortHandle = udp(Mstate.twoPhotonUDP,'RemotePort',7000);
% 
% set(TwoPcomState.serialPortHandle, 'OutputBufferSize', 1024)
% set(TwoPcomState.serialPortHandle, 'InputBufferSize', 1024)
% set(TwoPcomState.serialPortHandle, 'Datagramterminatemode', 'off')

%Establish serial port event callback criterion  
% TwoPcomState.serialPortHandle.BytesAvailableFcnMode = 'Terminator';
% TwoPcomState.serialPortHandle.Terminator = '~'; %Magic number to identify request from Stimulus ('c' as a string)
% TwoPcomState.serialPortHandle.bytesavailablefcn = @Displaycb;  

% open and check status 
fopen(TwoPcomState.serialPortHandle);
stat=get(TwoPcomState.serialPortHandle, 'Status');
if ~strcmp(stat, 'open')
    disp([' StimConfig: trouble opening port; cannot proceed']);
    TwoPcomState.serialPortHandle=[];
    out=1;
    return;
end
