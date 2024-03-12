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
        
     case 'VPX'
        
        % VPixx
        
        %Display size
        Mstate.screenXcm = 53.2;  %mm 
        Mstate.screenYcm = 30;  

        load('~/Desktop/calibration_stuff/measurements/Vpixx/LUT.mat','bufLUT')
        
        
        setImTrx([],[],[]);
        
        
    case 'LIN'   %load a linear table
        
        Mstate.screenXcm = 32.5;
        Mstate.screenYcm = 24;        
        
        bufLUT = (0:255)/255;
        bufLUT = bufLUT'*[1 1 1];
        
        setImTrx([],[],[]); 
        
   case 'Proj'   %load a linear table
        
        %Aspect ratio is 1.6 for DLP4500
%         Mstate.screenXcm = 25*2.54;
%         Mstate.screenYcm = 15.5*2.54;
         Mstate.screenXcm = 44; % 6/11/18 cropped size is 45.5
         Mstate.screenYcm = 28; % 6/11/18 cropped size is 27

        
        bufLUT = (0:255)/255;
        bufLUT = bufLUT'*[1 1 1];        
 
        
        %OutputPts are the location of the purple pts within the green grid
        %load('~/Desktop/calibration_stuff/GeometryCorrection/UVproj 6_7_18 close.mat','xyin','xyout') % proj screen is fully pushed in
        
        
        setImTrx([],[],[]);  %2 is for the "green gun", which is what the code transforms
        
        
	case 'Proj2'   %load a linear table
        
        %Aspect ratio is 1.6 for DLP4500; remember to update every time when you move
        %projector or screen locations
        Mstate.screenXcm = 44;
        Mstate.screenYcm = 28;
        
        bufLUT = (0:255)/255;
        bufLUT = bufLUT'*[1 1 1];        
 
        
        %OutputPts are the location of the purple pts within the green grid
        load('~/Desktop/calibration_stuff/GeometryCorrection/Dallas12_07_18','xyin','xyout')
            
         %% Fine adjustment of small shifts (units are pixels)
         Yshift = 5; %Positive moves green to the right
         Xshift = 0.5; %negative moves UV up

        xyout(:,2) = xyout(:,2)+Yshift;
        xyout(:,1) = xyout(:,1)+Xshift;
        
        setImTrx(xyin,xyout,2); %2 is for the "green gun", which is what the code transforms
        
        %setImTrx([],[],[]); 
        
       
end


Screen('LoadNormalizedGammaTable', screenPTR, bufLUT);  %gamma LUT

