function rgb = getRGBGain(theta,phi,Mag,contrastNorm)

%input is direction in S/M opsin space, but contrast (R) is not actual cone
%contrast that is normalized by the gray value

%If contrastNorm = 1, it normalizes to have Weber contrast dF/mean
%normalized.  I'm not sure why I would do it any other way, but I was.

%rgb(1) is the 'redgain'
%rgb(2) is the 'greengain'
%rgb(3) is the 'bluegain'

%Ti is the matrix that comes from inverting the dot product between the
%cone sensitivity functions and the monitor spectra;
%i.e. the inverse of: T = [Mcone(:) Scone(:)]'*[green(:) UV(:)]; 

Ti = [56.6927  -54.8512    6.0840;
   -9.0458   25.1595   -7.1038;
   -0.6170   -0.6725   18.2526] % 5/9/16

if contrastNorm
    T=inv(Ti); %Bring it back to the noninverted matrix    
    alpha = sum(T,2);
    Ti(:,1) = Ti(:,1)*alpha(1); %Correct for the different gray points of the cones
    Ti(:,2) = Ti(:,2)*alpha(2);
    Ti(:,3) = Ti(:,3)*alpha(3);
end

%Set the maximum amount of gain:
normer = sqrt(sum(Ti.^2,2));
normer = max(normer);
Ti(1,:) = Ti(1,:)/normer;
Ti(2,:) = Ti(2,:)/normer;
Ti(3,:) = Ti(3,:)/normer;

S = Mag*sin(phi); %S is z-axis
LMmag = Mag*cos(phi);
L = LMmag*cos(theta*pi/180); %L is x-axis
M = LMmag*sin(theta*pi/180); %M is y-axis

rgb = Ti*[L M S]';

