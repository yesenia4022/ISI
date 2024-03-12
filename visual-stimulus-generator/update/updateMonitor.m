function updateMonitor


global Mstate screenPTR 


switch Mstate.monitor
    
    case 'LCD'
        
        Mstate.screenXcm = 52.5;  
        Mstate.screenYcm = 29.5;       
      
        load('~/Desktop/calibration_stuff/measurements/LCD 2-15-16 PR655/LUT.mat','bufLUT')
        
        setImTrx([],[],[]);
        
        %bufLUT = (0:255)/255;
        %bufLUT = bufLUT'*[1 1 1];
        
    case 'CRT'
        
        % CRT; NEC 'BRIGHT SCREEN' at UT
        
        %Display size
        Mstate.screenXcm = 16*2.54;  %mm
        Mstate.screenYcm = 12*2.54;  

        load('~/Desktop/calibration_stuff/measurements/CRT NEC 8-28-15/LUT.mat','bufLUT')
        
        setImTrx([],[],[]);
        
    case 'LIN'   %load a linear table
        
        Mstate.screenXcm = 47.6;
        Mstate.screenYcm = 26.7;        
        
        bufLUT = (0:255)/255;
        bufLUT = bufLUT'*[1 1 1];
        
        setImTrx([],[],[]); 
        
    case 'Proj'   %load a linear table
        
        %Aspect ratio is 1.6 for DLP4500
%         Mstate.screenXcm = 25*2.54;
%         Mstate.screenYcm = 15.5*2.54;
         Mstate.screenXcm = 42.55; % 6/11/18 cropped size is 45.5
         Mstate.screenYcm = 26.67; % 6/11/18 cropped size is 27

        
        bufLUT = (0:255)/255;
        bufLUT = bufLUT'*[1 1 1];        
 
        
        %OutputPts are the location of the purple pts within the green grid
        load('~/Desktop/calibration_stuff/GeometryCorrection/NewProjHouse 10_25_19','xyin','xyout') % proj screen is fully pushed in
        
	case 'Proj2'   %load a linear table
            %Proj2 is for telescoped image- Gaby 10/25/19
        
        %Aspect ratio is 1.6 for DLP4500; remember to update every time when you move
        %projector or screen locations
         Mstate.screenXcm = 64  %5/17/20
%                             63.55; % 10/25/19- Gaby new size of telescoped is 63.55
                                    %6/11/18 cropped size is 45.5
         Mstate.screenYcm = 40   %5/17/20
%                              39.55; % 10/25/19- Gaby new size of telescoped is 39.55
                                    %6/11/18 cropped size is 27
        
        bufLUT = (0:255)/255;
        bufLUT = bufLUT'*[1 1 1];        
 
        setImTrx([],[],[]); 
        
        load('~/Desktop/calibration_stuff/GeometryCorrection/UVproj 01_10_20_full','xyin','xyout') % proj screen is fully pushed in

end


Screen('LoadNormalizedGammaTable', screenPTR, bufLUT);  %gamma LUT

