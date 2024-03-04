function LFPout = BandErase(LFP,fs,hl,hh)

%This one has the line noise cancelation as optional (lnflag).  Also, you
%can put in 0 for hh if you don't wan't a high pass filter, and/or inf for
%hl if you don't want a low-pass filter.

%This one takes in the channels (i.e. not organized by space) and uses a
%band-pass Butterworth... It also uses a notch filter.

%hh - high-pass cutoff in Hz
%hl - low-pass cutoff in Hz

N = length(LFP(1,:));

w = linspace(-pi,pi,N);

ordL = 8;
ordH = 8;
wl = 2*pi*hl/fs;
wh = 2*pi*hh/fs;
Butter_LP = 1./(1+(1i*w/wl).^(2*ordL));
Butter_HP = 1./(1+(wh./(1i*w)).^(2*ordH));


if hh ~= 0 & hl == inf
    Butter = 1-Butter_HP;
elseif hh == 0 & hl ~= inf
    Butter = 1-Butter_LP;
elseif hh ~= 0 & hl ~= inf
    Butter = 1-Butter_LP.*Butter_HP;
else 
    Butter = ones(size(Butter_LP));
end


H = Butter;



%%%Show Filter%%%
% co = -3;  
% id = find(10*log10(H) >= co);
% figure,plot(fs*w(id(1):id(end))/(pi*2),10*log10(H(id(1):id(end)))), ylabel('dB'), xlabel('Hz')
% ylim([co 1])
H = fftshift(H);  %make it go 0 to 2pi
%figure,plot(w*fs/(2*pi),fftshift(real(ifft(H))))

% figure,plot(w*fs/(2*pi),fftshift(H))
% hold on, plot(w*fs/(2*pi),fftshift(abs(fft(zscore(LFP)))),'r')
% asdf
%%%%%%%%%%%%%%%%%%%%

LFPout = NaN*zeros(size(LFP));

%figure,plot(fftshift(w)*fs/(2*pi),abs(fft(LFP(38,:))),'r')

for i = 1:length(LFP(:,1))
    s = LFP(i,:);
    s = real(ifft(abs(H).*fft(s)));
    LFPout(i,:) = s;
end


%hold on,plot(fftshift(w)*fs/(2*pi),abs(fft(LFPout(1,:))))
%asdf
