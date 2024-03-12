function updateMonitorValues

global Mstate

%Putting in the right pixels is not important because the stimulus computer
%asks for the actual value anyway.  It only matters if the analysis needs
%the right number of pixels (like retinotopy stimuli).

switch Mstate.monitor
    
    case 'LCD' 
        
        Mstate.screenXcm = 52.5;
        Mstate.screenYcm = 29.5;
        Mstate.xpixels = 1024;
        Mstate.ypixels = 768;
        
    case 'Proj'
        
        %Mstate.screenXcm = 46.7;
        %Mstate.screenYcm = 29.8;
        
%         Mstate.screenXcm = 42.55;
%         Mstate.screenYcm = 26.67;
        Mstate.screenXcm = 41.5; %7/23/18
        Mstate.screenYcm = 26.5; %7/23/18
        Mstate.xpixels = 912;
        Mstate.ypixels = 1140;
        
%         Mstate.T_MSGU = [0.1571 0.0695; 0.0002 0.0328];  %outer product of cones x gun spectra in dichromat: [Mcone Scone]'*[Green UV];
%         Mstate.CalibrationNotes = 'Calibrated 3/26/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated because projectors were moved by allison, unclear what the movement was';
%         Mstate.T_MSGU = [0.1837 0.0619; 0.0002 0.0287; 0.1607 0.0686];  %outer product of cones x gun spectra in dichromat: [Mcone Scone rod]'*[Green UV];
%         Mstate.CalibrationNotes = 'Calibrated 6/07/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated as a sanity check. Buffer values both at 128';
        Mstate.T_MSGU = [0.1823 0.0898; 0.0013 0.0587; 0.1600 0.0985];  %outer product of cones x gun spectra in dichromat: [Mcone Scone rod]'*[Green UV];
        Mstate.CalibrationNotes = 'Calibrated 7/23/18. Old EKB Projector, NEW 405 LED. New EKB fly eye, and condensor from Dallas projector. Replaced burnt fly eye from 7/06/18. 0.01" Teflon. Green 100% offset. Buffer values both at 128';

    case 'Proj2'
        
%         Mstate.screenXcm = 77; %previous screen size
%         Mstate.screenYcm = 48.5; %previous screen size
%         Mstate.screenXcm = 42.55;
%         Mstate.screenYcm = 26.67;
        Mstate.screenXcm = 41.5; %7/23/18
        Mstate.screenYcm = 26.5; %7/23/18
        Mstate.xpixels = 912;
        Mstate.ypixels = 1140;
       
        %These calibration values are not correct.  I just copied them for
        %the closer projector position above
%        Mstate.T_MSGU = [0.1571 0.0695; 0.0002 0.0328];  %outer product of cones x gun spectra in dichromat: [Mcone Scone]'*[Green UV];
%        Mstate.CalibrationNotes = 'Calibrated 3/26/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated because projectors were moved by allison, unclear what the movement was';
        Mstate.T_MSGU = [0.1823 0.0898; 0.0013 0.0587; 0.1600 0.0985];  %outer product of cones x gun spectra in dichromat: [Mcone Scone rod]'*[Green UV];
        Mstate.CalibrationNotes = 'Calibrated 7/23/18. Old EKB Projector, NEW 405 LED. New EKB fly eye, and condensor from Dallas projector. Replaced burnt fly eye from 7/06/18. 0.01" Teflon. Green 100% offset. Buffer values both at 128';

    case 'CRT'

%         Mstate.screenXcm = 32.5;
%         Mstate.screenYcm = 24;
        
        Mstate.screenXcm = 16*2.54;
        Mstate.screenYcm = 12*2.54;
        
        Mstate.xpixels = 1024;
        Mstate.ypixels = 768;
        
     case 'TEL'

%         Mstate.screenXcm = 32.5;
%         Mstate.screenYcm = 24;
        
        Mstate.screenXcm = 121;
        Mstate.screenYcm = 68.3;
        
        Mstate.xpixels = 1024;
        Mstate.ypixels = 768;
        
end