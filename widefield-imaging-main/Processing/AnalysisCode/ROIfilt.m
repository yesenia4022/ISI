function imout = ROIfilt(im,bw,LPsig,HPdiameter, varargin)

M = size(im,1);
N = size(im,2);

%Smooth (low pass)
if LPsig
    hh = fspecial('gaussian',[M N],LPsig);
    bwSmear = ifft2(fft2(bw).*abs(fft2(hh)));
    bwSmear = bwSmear .* bw;
    im = ifft2(fft2(im.*bw).*abs(fft2(hh))) ./ bwSmear;  %Only smear within the ROI, 'bw'
    im(~bw) = 0;
end

%Subtract local mean (high pass)
if ~isinf(HPdiameter)
    [x y] = meshgrid(1:N,1:M);
    r = sqrt((x-mean(x(:))).^2  + (y-mean(y(:))).^2);
    hh = zeros(M,N);
    hh(find(r<=HPdiameter/2)) = 1; hh = hh/sum(hh(:));  %Make local mean filter
    
    if ~isempty(varargin{1})
        hh = fspecial('gaussian',[M N],HPdiameter);
    end
    
    bwSmear = ifft2(fft2(bw).*abs(fft2(hh)));
    bwSmear = bwSmear .* bw;
    imdum = ifft2(fft2(im.*bw).*abs(fft2(hh))) ./ bwSmear;  %Only smear within the ROI, bw
    imdum(~bw) = 0;
    imout = im-imdum;
else
    imout = im;
end