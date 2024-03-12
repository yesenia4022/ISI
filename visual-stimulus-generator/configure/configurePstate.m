function configurePstate(modID)


switch modID
    
    case 'PG'
        
        configPstate_perGrater
        
    case 'FG'
      
        configPstate_flashGrater
        
    case 'RD'
        
        configPstate_Rain
        
    case 'FN'
        
        configPstate_Noise
        
    case 'MP'
        
        configPstate_Mapper
        
     case 'CM'
        
        configPstate_cohMotion
        
     case 'GC'
        
        configPstate_geoCalibration
        
     case 'SI'
        
        configPstate_savedImages

     case 'NM'
        
        configPstate_naturalMovie
        
     case 'US'
        
        configPstate_uncertaintyStimulus
        
end

