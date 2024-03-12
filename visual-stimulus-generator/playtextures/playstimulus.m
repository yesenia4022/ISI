function playstimulus(modID)


HideCursor

switch modID
    
    case 'PG'  %Periodic Grater        

        if getParamVal('mouse_bit');
            playmanualgrater
        else
            playgrating_periodic
        end
                
    case 'FG'  %Flash Grater
        
        if getParamVal('FourierBit')
            playgrating_flashHartley
        else
            playgrating_flashCartesian
        end
        
    case 'RD'  %Raindropper
        
        playrain
        
    case 'FN'  %Filtered Noise
        
        playnoise
        
    case 'MP'  %Mapper
        
        playmapper
        
    case 'CM'  %Coherent Motion
        
        playcohmotion
        
    case 'GC' %Geometric calibration
        
        GeoCalibration
        
    case 'SI' %Saved images 
        
        playSavedImages
        
    case 'NM' % Natural Movie
        
        playnaturalmovie
        
    case 'US' % yb: uncertainty stimulus
        
        playUncertaintyStimulus
        
end

ShowCursor
    