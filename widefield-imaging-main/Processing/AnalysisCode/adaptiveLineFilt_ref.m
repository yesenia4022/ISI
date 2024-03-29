function [sigfilt] = adaptiveLineFilt_ref(signal,ref1,ref2,f0,sp,W,shiftN,varargin)

%_ref uses a reference line to identify noise

if isempty(varargin)
    [e,v] = dpss(W,1.5);
else
    e = varargin{1};
end

%lazy way of computing the number of shifts
N = length(signal);
for i = 1:N

    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;

    if stdom(end)>N
        break
    end

    Nshift = i;

end

freqR = [f0-.1*f0  f0+.1*f0];  %Allowed drift in the peak


%Generate time frequency matrices

win = hann(W);
%win = ones(W,1);

lineN = zeros(size(signal));
winMod = zeros(size(signal));

fs = 1000/sp;
for i = 1:Nshift
    
    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;
    
    n = 0:(W-1);

    sigdum = signal(stdom);
    sigdum = sigdum-mean(sigdum);
    sigdum = sigdum.*win;
    
    ref1dum = ref1(stdom);
    ref1dum = ref1dum-mean(ref1dum);
    ref1dum = ref1dum.*win;
    ref2dum = ref2(stdom);
    ref2dum = ref2dum-mean(ref2dum);
    ref2dum = ref2dum.*win;
    
    [ref1pow fdom] = getPowerSpect(ref1dum,sp,W,round(W/2));
    [ref2pow fdom] = getPowerSpect(ref2dum,sp,W,round(W/2));
    [sigpow fdom] = getPowerSpect(sigdum,sp,W,round(W/2));
    
    Xpow = abs(ref1pow.*ref2pow.*sigpow);
    
    [dum id1] = min(abs(fdom-freqR(1)));
    [dum id2] = min(abs(fdom-freqR(2)));
    [dum id] = max(Xpow(id1:id2));
    id = id+id1-1;
    f1(i) = fdom(id);
    
    F1 = sigdum(:)'*exp(1i*2*pi*f1(i)*n(:)/fs);

    tcF1_line = abs(F1)*cos(2*pi*f1(i)*n(:)/fs - angle(F1));

    amp = sum((cos(2*pi*f1(i)*n(:)/fs).^2).*win);
    tcF1_line = tcF1_line/amp;
    tcF1_line = tcF1_line.*win;
    
    lineN(stdom) = lineN(stdom) + tcF1_line;
    
    winMod(stdom) = winMod(stdom) + win;

end
id = find(winMod == 0);
winMod(id) = 1;
lineN = lineN./winMod;

sigfilt = signal-lineN;
