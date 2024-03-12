function configureMstate

global Mstate

NHPflag =0;

Mstate.anim = 'xx0';
Mstate.unit = '000';
Mstate.expt = '000';

Mstate.hemi = 'right';

Mstate.running = 0;

if NHPflag
    
    Mstate.screenDist = 60;
    
    Mstate.monitor = 'CRT';  %This should match the default at the master. Otherwise, they will differ, but only at startup
    
    %'updateMonitor.m' happens in 'screenconfig.m' at startup
    
    Mstate.syncSize = 4;  %cm
    
    
else
    
    Mstate.screenDist = 10;
    
    Mstate.monitor = 'LIN';  %This should match the default at the master. Otherwise, they will differ, but only at startup
    
    %'updateMonitor.m' happens in 'screenconfig.m' at startup
    
    Mstate.syncSize = .5;  %cm
    
end


