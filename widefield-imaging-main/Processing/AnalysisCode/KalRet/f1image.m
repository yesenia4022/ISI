function y = f1image(cond)

% Compute the F1 for one condition

global ACQinfo

acqPeriod = ACQinfo.linesPerFrame*ACQinfo.msPerLine/1000;  %sec per acquired frame
Npix = ACQinfo.linesPerFrame*ACQinfo.pixelsPerLine;
Fs = 1000*ACQinfo.pixelsPerLine/ACQinfo.msPerLine;

nr = getnorepeats(1);
%for r = 0:0
for r = 1:nr
    r
    clear CHs
    CHsSync = GetTrialData([0 0 1 0],[cond r]);
    Nframes = length(CHsSync{1}(1,1,:));
    
    dim = size(CHsSync{1});
    synccourse = NaN*ones(prod(dim),1);
    for i = 1:Nframes
        dum = CHsSync{1}(:,:,i)';
        synccourse(1+(i-1)*dim(2)*dim(1):i*dim(2)*dim(1)) = dum(:);
    end
    
    clear CHsSync
    
    synctimes = getLCDsynctimes(synccourse,Fs);
    synctimes = synctimes(2:end-1);  %First and last synces are before and after pre/post delay
    
    clear synccourse dum
    
    T = mean(diff(synctimes(1:end-1)));  %Don't use last sync here because the last cycle may be incomplete 
    
    framest = 0:Nframes-1;
    framest = framest*acqPeriod;
    frameang = NaN*ones(1,length(framest));
    
    for i = 1:length(synctimes)-1  %The last sync signifies stimulus end
        id = find(framest>=synctimes(i) & framest<synctimes(i)+T);
        frameang(id) = 2*pi*(framest(id) - synctimes(i))/T;
    end

    idx = find(~isnan(frameang));
    
    CHs = GetTrialData([1 0 0 0],[cond r]);
    
    muIm = mean(CHs{1}(:,:,idx),3);

    tic
    for j=1:length(idx)

        img = CHs{1}(:,:,idx(j)) - muIm;  %Subtract f0 leakage
     
        if j==1
            acc = zeros(size(img));
        end
        
        %Don't put a negative sign in the exponent if you want higher
        %phases to mean a rightward shift (i.e. greater delay) of the waveform.
        %Apparently, the Fourier trx assumes the opposite
        acc = acc + exp(1i*frameang(idx(j))).*img; 
        
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

   y{r} = acc;

end

    
    