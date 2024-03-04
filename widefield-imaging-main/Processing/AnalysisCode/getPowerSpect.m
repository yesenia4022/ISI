function [sigpow f_dom] = getPowerSpect(sig,sp,W,shiftN,varargin)

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

HH = zeros(W,shiftN*length(e(1,:))); %preallocate
z = 1;
for i = 1:Nshift

    stdom = 1+(i-1)*shiftN:(i-1)*shiftN+W;
    
    sigdum = sig(stdom);
    sigdum = sigdum-mean(sigdum);

    for s = 1:length(e(1,:))
        HH(:,z) = fft(win.*sigdum.*e(:,s));

        z = z+1;
    end

end


sigpow = sum(abs(HH).^2,2); %power spectrum

f_dom = linspace(0,1000/sp,length(sigpow)+1);
f_dom = f_dom(1:end-1);

sigpow = sigpow(1:round(length(sigpow)/2));
f_dom = f_dom(1:round(length(f_dom)/2));


% figure,plot(f_dom,angle(coh))
% figure,plot(f_dom,abs(coh))