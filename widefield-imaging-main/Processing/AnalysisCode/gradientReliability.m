function GradSel = gradientReliability(im,h)

[dIdy dIdx] = gradient(im);

GradSel = sqrt(dIdx.^2 + dIdy.^2);


gradvec = (dIdx + 1i*dIdy)./sqrt(dIdx.^2 + dIdy.^2);  %Vector image of unity amplitude
%gradvec = (dIdx + 1i*dIdy);
%gradmag = sqrt(dIdx.^2 + dIdy.^2);

id = find(isnan(gradvec)); gradvec(id) = 0;
%id = find(isnan(gradmag)); gradmag(id) = 0;


h = h/sum(h(:));
GradSel = abs(ifft2(fft2(gradvec).*abs(fft2(h))));  %Gradient selectivity index (0 to 1)

% GradMag = abs(ifft2(fft2(gradmag).*abs(fft2(h))));  %Gradient selectivity index (0 to 1)
% GradSel = GradSel./GradMag;