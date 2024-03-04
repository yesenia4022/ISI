function yF = getResPeaks(y,sp,winT)

fdom = linspace(0,1/(sp/1000),length(y)+1);
fdom = fdom(1:end-1);

sig = round(winT/sp/2);
Wshift = sig;
idom = 0:length(y)-1;
yF = 0;
for i = 1:round(length(y)/Wshift)
    
    shift = (i-1)*Wshift;
    Win = exp(-(idom-shift).^2/(2*sig^2));
    yF = yF + abs(fft(Win.*y));
    
end

%figure,plot(fdom,yF)




