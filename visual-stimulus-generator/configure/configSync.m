function configSync

global daq

%daq = DaqDeviceIndex; %This does not work on new OS

daq = 4; %A guess  Do D = PsychHID('Devices'); D.product; to find an index
if ~isempty(daq)
    
    DaqDConfigPort(daq,0,0);    
    
    %DaqDOut(daq, 0, 0); 
    
else
    
    'Daq device does not appear to be connected'
    
end