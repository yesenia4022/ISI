function [Xpow sig1pow sig2pow f_dom] = getXPowerSpect(sig1,sig2,sp,W,shiftN,varargin)

if isempty(varargin)
    [e,v] = dpss(W,1.5);
else
    e = varargin{1};
end

%Generate time frequency matrices

%lazy way of computing the number of shifts
N = length(sig);
for i = 1:N

    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;

    if stdom(end)>N
        break
    end

    Nshift = i;

end

win = hann(W);
%win = ones(W,1);

HH1 = zeros(W,shiftN*length(e(1,:))); %preallocate
HH2 = zeros(W,shiftN*length(e(1,:))); %preallocate

z = 1;
for i = 1:Nshift

    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;
    
    sig1dum = sig1(stdom);
    sig1dum = sig1dum-mean(sig1gdum);
    sig2dum = sig2(stdom);
    sig2dum = sig2dum-mean(sig2dum);

    for s = 1:length(e(1,:))
        HH1(:,z) = fft(win.*sig1dum.*e(:,s));
        HH2(:,z) = fft(win.*sig2dum.*e(:,s));

        z = z+1;
    end

end


sig1pow = mean(abs(HH1).^2,2); %power spectrum
sig2pow = mean(abs(HH2).^2,2); %power spectrum

Xpow = mean(HH1.*conj(HH2),2);

f_dom = linspace(0,1000/sp,length(sig)+1);
f_dom = f_dom(1:end-1);

sig1pow = sigpow(1:round(length(sigpow)/2));
sig2pow = sigpow(1:round(length(sigpow)/2));
Xpow = sigpow(1:round(length(sigpow)/2));
f_dom = f_dom(1:round(length(f_dom)/2));


% figure,plot(f_dom,angle(coh))
% figure,plot(f_dom,abs(coh))