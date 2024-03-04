function synctimes = getFlashGraterSynctimes_CRT(syncwave,Fs)

%Produces a vector corresponding to the rising edge times of syncwave

%Decimate, because otherwise it sometimes runs out of memory at the fft line
Dec = 40;  
syncwave = syncwave(Dec:Dec:end);
Fs = Fs/Dec;

%%%%%%%%%%%%%%%%
%We are getting annoying fucking zeros during the bidir scan.  This gets rid of them: 
id = find(syncwave == 0);
for i = 1:length(id)    
    id1 = max(id(i)-20,1);
    id2 = min(id(i)+20,length(syncwave));
    syncwave(id(i)) = mean([syncwave(id1) syncwave(id2)]);
end
syncwave = medfilt1(syncwave,3); %they are not all zero this gets rid of the rest
%%%%%%%%%%%%%%%%%%%%%%%%%%%

syncwaveF = fft(syncwave-mean(syncwave));

dF = Fs/length(syncwave);
milookup = 80;
idmi = round(milookup/dF);
[dum idma] = max(syncwaveF(idmi:end/2));
idma = idma + idmi-1; 

noiseF0 = (idma-1)*dF;
%noiseF0 = fdom(id);
%noiseF0 = 10;
W = round(Fs/noiseF0);

H = zeros(1,length(syncwave));
H(1:W) = 1/W;
H = abs(fft(H'));
%H = H.^2;

syncwave = ifft(fft(syncwave).*H);
high = max(syncwave);
low = min(syncwave);
thresh = (high+low)/2;

%%%
%thresh = 0.2;
%%%
syncwave = sign(syncwave-thresh);
id = find(syncwave == 0);
syncwave(id) = 1;
syncwave = abs(diff((syncwave+1)/2));

synctimes = find(syncwave == 1) + 1;
synctimes = synctimes/Fs;

clear syncwave

