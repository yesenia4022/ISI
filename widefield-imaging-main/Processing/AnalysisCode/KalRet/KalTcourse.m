function y = KalTcourse(cond)

%% Compute the F1 for one condition

global pepANA ACQinfo

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine/1000;  %sec per acquired frame
Npix = ACQinfo.linesPerFrame*ACQinfo.pixelsPerLine;

pepsetcondition(cond);
nr = pepgetnorepeats;
%for r = 0:0
for(r=0:nr-1)
    pepsetrepeat(r)
    r
    clear CHs
    CHsSync = GetTrialData([0 0 1]);
    Nframes = length(CHsSync{1}(1,1,:));
    
    dim = size(CHsSync{1});
    synccourse = NaN*ones(prod(dim),1);
    for i = 1:Nframes
        dum = CHsSync{1}(:,:,i)';
        synccourse(1+(i-1)*dim(2)*dim(1):i*dim(2)*dim(1)) = dum(:);
    end
    
    clear CHsSync
    
    synctimes = getsynctimes(synccourse)/1000;
    
    clear synccourse dum
    
    T = mean(diff(synctimes));
    
    framest = 0:Nframes-1;
    framest = framest*acqPeriod;
    frameang = NaN*ones(1,length(framest));
    
    for i = 1:length(synctimes)
        id = find(framest>=synctimes(i) & framest<synctimes(i)+T);
        frameang(id) = 2*pi*(framest(id) - synctimes(i))/T;
    end

    idx = find(~isnan(frameang));
    
    CHs = GetTrialData([1 0 0]);
    
    phasedom = 0:5:360-5;

    Tens = zeros(dim(1),dim(2),length(phasedom));
    counter = zeros(1,length(phasedom));
    tic
    for j=1:length(idx)
        
        if j ~= 1
            if frameang(idx(j)) < frameang(idx(j-1))
                iniImage = CHs{1}(:,:,idx(j));
            end
        else
            iniImage = CHs{1}(:,:,idx(j));
        end
        
        ang = frameang(idx(j))*180/pi;
        [dum iddom] = min(abs(phasedom-ang));
        
        img = CHs{1}(:,:,idx(j));
        
        Tens(:,:,iddom) = Tens(:,:,iddom) + img - iniImage;
        
        counter(iddom) = counter(iddom) + 1;
        
    end
    toc
    
    for j = 1:length(phasedom)
        Tens(:,:,j) = Tens(:,:,j)/counter(j);
    end
    
    y{r+1} = Tens;

end

    
    