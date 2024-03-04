function resetImage(trx)

global imstate

switch trx
    
    case 'UD'        
        imstate.imanat = flipud(imstate.imanat);    
        imstate.imanatFunc = flipud(imstate.imanatFunc);   
        imstate.imfunc = flipud(imstate.imfunc);  
        imstate.bw = flipud(imstate.bw);
        imstate.mag = flipud(imstate.mag);
        
        imstate.fmaps{1} = flipud(imstate.fmaps{1});
        imstate.fmaps{2} = flipud(imstate.fmaps{2});
        imstate.sigMag = flipud(imstate.sigMag);
        
        imstate.areaBounds = flipud(imstate.areaBounds);
    case 'LR'        
        imstate.imanat = fliplr(imstate.imanat); 
        imstate.imanatFunc = fliplr(imstate.imanatFunc);   
        imstate.imfunc = fliplr(imstate.imfunc); 
        imstate.bw = fliplr(imstate.bw);
        imstate.mag = fliplr(imstate.mag);
        
        imstate.fmaps{1} = fliplr(imstate.fmaps{1});
        imstate.fmaps{2} = fliplr(imstate.fmaps{2});
        imstate.sigMag = fliplr(imstate.sigMag);
        
        imstate.areaBounds = fliplr(imstate.areaBounds);
    case 'rotate'        
        imstate.imanat = flipud(imstate.imanat');  %Rotate 90deg cc
        imstate.imanatFunc = flipud(imstate.imanatFunc');  
        imstate.imfunc = flipud(imstate.imfunc');
        imstate.bw = flipud(imstate.bw');
        imstate.mag = flipud(imstate.mag');
        
        imstate.fmaps{1} = flipud(imstate.fmaps{1}');
        imstate.fmaps{2} = flipud(imstate.fmaps{2}');
        imstate.sigMag = flipud(imstate.sigMag');
        
        imstate.areaBounds = flipud(imstate.areaBounds');
end






