function configureMstate


NHPflag = 0;

global Mstate

Mstate.anim = 'xx0';
Mstate.unit = '000';
Mstate.expt = '000';

Mstate.hemi = 'left';

if NHPflag

    Mstate.screenDist = 60;
    Mstate.monitor = 'CRT';  %This should match the default value in Display
    Mstate.syncSize = 4;  %Size of the screen sync in cm
    
else
    
    Mstate.screenDist = 10;
    Mstate.monitor = 'LIN';  %This should match the default value in Display
    Mstate.syncSize = 2;  %Size of the screen sync in cm
    
end

updateMonitorValues

Mstate.running = 0;

Mstate.analyzerRoot = 'C:\AnalyzerFiles'; 

Mstate.analyzerRoot_2p = '\\2PACQUISITION\AnalyzerFiles';

Mstate.dataRoot = 'D:\ISIdata'; 

Mstate.stimulusUDP = '128.114.78.199';   %220 %Neighbor (ISI computer) Slave (Mac)

Mstate.twoPhotonUDP = '128.83.22.238';  %Neighbor (two-photon computer)

Mstate.twoP = 0;

Mstate.WF = 0;

Mstate.Ephys = 0;

Mstate.framerate = 60;

Mstate.electrode = 'A8x8';

Mstate.camera = 'Dalsa';

Mstate.framerate
