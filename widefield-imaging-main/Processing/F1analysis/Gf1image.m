function [y s] = Gf1image(cond,varargin)

%%% Compute the F1 for one condition

global Analyzer GF1_handles datadir

%Get acquired bit depth by pulling in a raw frame
ue = datadir(end-8:end-1);
fname = [datadir ue  '_' sprintf('%03d',0)];
fnamedum = [fname '_f1'];
load(fnamedum)
k = whos('im');
bitDepth = str2num(k.class(5:end));


inputs = varargin{1};

if length(inputs)==2    %if we want to analyze a pixel.
    N = length(inputs{1}(:,1));   %# of pixels;
    ps = inputs{2};  %pixel size
    pixels = inputs{1};  %Matrix of pixel locations; 1st column is y values
    for i = 1:N
        yr(i,:) = ((pixels(i,1):pixels(i,1)+ps-1))-floor(ps/2);   %yrange of first pixel
        xr(i,:) = ((pixels(i,2):pixels(i,2)+ps-1))-floor(ps/2);   %xrange of first pixel
    end
else
    s = [];
end

repSet = get(GF1_handles.reps,'string');
if strcmp(repSet,'All')
    repdom = 1:getnorepeats(1);
else
    repdom = eval(repSet);
end

frate = Analyzer.syncInfo{1}.frameRate;

for rid = 1:length(repdom)
    r = repdom(rid)
    
    tno = Analyzer.loops.conds{cond}.repeats{r}.trialno
    Analyzer.syncInfo{tno}.frameRate = 60
    Grabtimes = Analyzer.syncInfo{tno}.acqSyncs;
    Grabtimes = Grabtimes - Grabtimes(1);
    
    %% Old code
    %Stimulus starts on 2nd sync, and ends on the second to last.  I also
    %get rid of the last bar rotation (dispSyncs(end-1)) in case it is not an integer multiple
    %of the stimulus trial length
    %Disptimes = Analyzer.syncInfo{tno}.dispSyncs(2:end-1);
    %%
    
    preScreenFlips = round(getParamVal('predelay',0)*frate);
    preDelayEnd = preScreenFlips/frate;
    
    try
        periodTime = getParamVal('t_period',0)/frate;
    catch
        periodTime = getparam('block_Nimage')*getparam('h_per')/60*2;
    end
    Nperiods = floor(getParamVal('stim_time',0)/periodTime);
    
    Disptimes_start = [preDelayEnd + (0:Nperiods-1)*periodTime];
    Disptimes_end = [preDelayEnd + (1:Nperiods)*periodTime];
    %Disptimes = Disptimes(1:end-1); %throw out last one in case its not an integer multiple
    
    
    %T = getParamVal('t_period',0)/60;
    T = mean(diff(Disptimes_start)); %This one might be more accurate
    
    fidx = find(Grabtimes>Disptimes_start(1) & Grabtimes<Disptimes_end(end));  %frames during stimulus

    framest = Grabtimes(fidx)-Disptimes_start(1);  % frame sampling times in sec
    frameang = framest/T*2*pi;       %% convert to radians
    
    %Stack = loadTrialData(cond,r);
    
    k = 1;
   
    win = hann(length(fidx(1):fidx(end)));
    tic
    for j=fidx(1):fidx(end)
        
        if get(GF1_handles.negSignal,'value')
            img = 2^bitDepth-getTrialFrame(j,cond,r);
        else
            img = getTrialFrame(j,cond,r);
        end
        %img = 4096-double(Stack(:,:,j));
        
        if length(inputs)==2
            for i = 1:N
                sig(i,k) = mean2(img(yr(i,:),xr(i,:)));
            end
        end
     
        if j==fidx(1)
            acc = zeros(size(img));
            f0 =  zeros(size(img));
            ms =  zeros(size(img));
        end

        acc = acc + exp(1i*frameang(k)).*img;
        f0 = f0 + img;
        ms = ms + img.^2;

        k = k+1;

    end
    toc
    
  %% f0 = f0./(fidx(end)-fidx(1)+1);
   f0 = f0./(k-1); %mean
   ms = (ms./(k-1)); %mean-squared
   stdIm = (ms - f0.^2);
   acc = acc - f0*sum(exp(1i*frameang)); %Subtract f0 leakage
  %% acc = 2*acc ./ (length(fidx(1):fidx(end))); %Normalize for f1 amplitude
   acc = 2*acc ./ (k-1); %Normalize for f1 amplitude 
  
   %y{rid} = acc./f0;  %dF/F
   y{rid} = acc./stdIm; 
   
   if length(inputs)==2
       mp = mean(sig');             %find mean of each pixel in the repeat
       mp = meshgrid(mp,1:length(sig(1,:)))';
       s{rid} = sig;
       %s{rid} = sig-mp;
   end

end

    
    