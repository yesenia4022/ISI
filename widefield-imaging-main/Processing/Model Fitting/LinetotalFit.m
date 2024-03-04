function [thetabest bbest] = LinetotalFit(xx,yy,thetadom,bdom)

yy = yy(:);
xx = xx(:);


for thetaid = 1:length(thetadom)

    thet = thetadom(thetaid);
    m = tan(thet*pi/180);

    for bid = 1:length(bdom)

        b = bdom(bid);

        if thet < 45            
            dy = yy - (m*xx + b);
            dy = abs(dy);
            E = dy*cos(thet*pi/180);
        else
            dx = xx - (yy-b)/m;
            dx = abs(dx);
            E = dx*sin(thet*pi/180);
        end

        err(thetaid,bid) = sum(abs(E));

    end
end

smoother = zeros(size(err));
smoother(1:30) = hann(30);
err = squeeze(err);
err  = ifft(abs(fft(smoother)).*fft(err));


[thetabest bbest] = find(err == min(err(:)));
thetabest = thetadom(thetabest);
bbest = bdom(bbest);