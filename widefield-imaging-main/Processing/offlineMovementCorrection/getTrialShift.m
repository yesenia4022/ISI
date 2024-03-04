function [Vx Vy] = getTrialShift(im,temp)

%Light smoothing of the data across time
% tkern = [0 1 0];
% tkern = tkern/sum(tkern);
% kern = zeros(Idim(1),Idim(2),length(tkern));
% for i = 1:length(tkern)
%     kern(:,:,i) = ones(Idim(1),Idim(2))*tkern(i);
% end
% 
% smoother = zeros(size(CH));
% smoother(:,:,1:length(tkern)) = kern;
% 
% smoother = abs(fft(smoother,[],3));
% CHsfilt = ifft(fft(CH,[],3).*smoother,[],3);

CH{1} = temp;
CH{2} = im;

%Smooth data in space
skern = hann(5)';
skern = skern'*skern;
skern = skern/sum(skern(:));


for k = 1:2
    CH{k}(1,:,:) = CH{k}(2,:,:);
    CH{k}(end,:,:) = CH{k}(end-1,:,:);
    CH{k}(:,1,:) = CH{k}(:,2,:);
    CH{k}(:,end,:) = CH{k}(:,end-2,:);
    CH{k}(:,end-1,:) = CH{k}(:,end-2,:);
    %CHsfilt(end-2:end,:) = ones(3,1)*CHsfilt(end-3,:);
    
    bw = ones(size(CH{k}));
    id = find(isnan(CH{k}));
    bw(id) = 0;
    
    CH{k} = roifilt2(skern,CH{k},bw);
    %CH{k} = ifft2(fft2(CH{k}).*smoother);

end



[dFdx1 dFdy1] = gradient(CH{2}); %change per pixel, and change per frame
[dFdx2 dFdy2] = gradient(CH{1}); %change per pixel, and change per frame

dFdx = (dFdx1+dFdx2)/2;
dFdy = (dFdy1+dFdy2)/2;

dFdt = CH{2} - CH{1};

xdum = dFdx(2:end-1,2:end-2);
ydum = dFdy(2:end-1,2:end-2);
tdum = dFdt(2:end-1,2:end-2);


dum = xdum.*ydum.*tdum;
id = find(~isnan(dum(:)));
pc = princomp([xdum(id) ydum(id) tdum(id)]);

H = [xdum(id) ydum(id)];  %plane should go through the origin: time derivative is zero when spatial derivs are zero
Vxy = inv(H'*H)*H'*tdum(id);

% pc12 = pc(:,1:2)';
% H = pc12(:,1:2); y = pc12(:,3);
% Vxy = inv(H'*H)*H'*y;

Vx = Vxy(1);
Vy = Vxy(2);

%figure,scatter(ydum(id),tdum(id),'.')
