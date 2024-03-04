function y = f1imageIT(cond)

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
    
    synccourse = [];
    for i = 1:Nframes
        dum = CHsSync{1}(:,:,i)';
        synccourse = [synccourse; dum(:)];
    end
    
    synctimes = getsynctimes(synccourse)/1000;

    clear CHsSync synccourse dum
    
    T = mean(diff(synctimes));
    
    framest = 0:Nframes-1;
    framest = framest*acqPeriod;
    frameang = NaN*ones(1,length(framest));
    
    for i = 1:length(synctimes)
        id = find(framest>=synctimes(i) & framest<synctimes(i)+T);
        frameang(id) = 2*pi*(framest(id) - synctimes(i))/T;
    end

%     frames = pepgetframeidx(ttag,[0 sync(end)]);
%     framest = pepANA.imaging.isync(frames(1):frames(2));  %integer sample values
%     framest = (framest - ts(1))/30;  %% frame sampling times in msec  (Cerebrus: 30 samp/ms)
%     frameang = framest/T*2*pi;       %% convert to radians

    idx = find(~isnan(frameang));
    
    CHs = GetTrialData([1 0 0]);
    
    muIm = mean(CHs{1}(:,:,idx),3);

    tic
    for j=1:length(idx)

        img = CHs{1}(:,:,idx(j)) - muIm;  %Subtract f0 leakage
     
        if j==1
            acc = zeros(size(img));
        end

        acc = acc + exp(-1i*frameang(idx(j))).*img;
        
        %accbias = accbias + exp(-1i*frameang(idx(j)));

    end
    toc
    
    pixperimage = ACQinfo.linesPerFrame*ACQinfo.pixelsPerLine;
    secperpix = ACQinfo.msPerLine/1000/ACQinfo.pixelsPerLine;
    delaymask = 0:pixperimage-1;
    delaymask = reshape(delaymask',ACQinfo.pixelsPerLine,ACQinfo.linesPerFrame)';
    delaymask = delaymask*secperpix/T*2*pi;  %Scan delay for each pixel in radians

    acc = acc.*exp(1i*delaymask);

  %% acc = 2*acc ./ (length(frames(1):frames(2))); %Normalize for f1 amplitude
   acc = 2*acc ./ length(idx); %Normalize for f1 amplitude 

   y{r+1} = acc;

end

    
    