function moveShutter(eye,pos)

global daq shutterState

LRports = [3 2];  %first element is lefteye port on 1208FS, second element is righteye port [0 to 7]

%eye will be a 0 or 1, so LRports(eye+1) will id the hardware port

shutterState=bitset(shutterState,LRports(eye),1-pos); %set 'eye' bit 0 or 1, to pos 0 or 1

disp(shutterState)

%DaqDOut(daq, 1, shutterState);