function gb = getGreenBlueGain(theta,R,contrastNorm)

%input is direction in S/M opsin space, but contrast (R) is not actual cone
%contrast that is normalized by the gray value

%If contrastNorm = 1, it normalizes to have Weber contrast dF/mean
%normalized.  I'm not sure why I would do it any other way, but I was.

%gb(1) is the 'greengain'
%gb(2) is the 'bluegain'

%"Tmatrix" comes from inverting the dot product between the
%cone sensitivity functions and the monitor spectra;
%i.e. the inverse of: T = [Mcone(:) Scone(:)]'*[green(:) UV(:)]; 


%Tmatrix = [12.7805 -23.5983; -0.0118 25.1847]; %Calibrated 4/20/16


%Tmatrix = [21.1055 -38.8909; -0.0297 25.6989]; %'Middle' Calibrated 5/8/16 
%Tmatrix = [35.6516 -65.9533; -0.1890 38.9525]; %'Right' Calibrated 5/8/16


%Tmatrix = [19.0449  -34.7339; -0.1270   32.1627]; %Calibrated 12/12/16 

%Tmatrix = [17.3483  -32.6743; -0.1070   61.5954]; %Calibrated 2/18/17 (New RGB, same UV, bigger display)

% Tmatrix = [18.2849  -34.0223; -0.1570   69.3223]; %Calibrated 3/29/17 (We noticed drop in UV power)

% Tmatrix = [10.5291  -28.6346; -0.0082   50.2133]; %Calibrated 2/24/17 "New projectors" . Has Keynote nUV and 100% offset RGB
%Tmatrix = [11.3075  -30.3970; -0.1856   58.2362]; %Calibrated 4/11/17 Keynote projectors from 1.164 while previous UV proj. getting repaired
% Tmatrix = [45.2300 -121.5880; -0.1856   58.2362];

%Tmatrix = [9.2651  -28.5003; -0.1912  54.0682]; %Calibrated 4/17/17. Keynote UV and 100% offset RGB.  UV current at 250. Green current at 160. 


%Tmatrix =  [15.5112  -6.9849; -0.0311  9.0630]; %Calibrated 6/20/17. Moved projectors closer .01" Teflon. New UV projector (EKB) with 385 LED. 0% offset RGB.  Green at 64, UV at 128
%Tmatrix = [22.3927 -10.7496; -.0082 19.7715]; % %Calibrated 6/23/17.  Dallas projector (385 LED) and 100% offset RGB in 1.164. Teflon.  current at 100 and 64 for UV and RGB, respectively
%Tmatrix = [25.8790  -56.9882;   -0.2355   47.3623]; % % Calibrated 7/14/17. Old EKB Projector, 405 LED. UV current 128, Green Current 64. 0.01" Teflon. Green 100% offset
%Tmatrix = [13.8724  -30.1537;   -0.0454   25.1747];%% Calibrated 10/16/17. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Projector was pushed up by 21.8cm new screen size is 46.7 by 29.8 cm
%Tmatrix = [6.7382  -14.2522; -0.0305   25.4097];%% Calibrated 3/02/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated due to sanity check of projector output, dicovered that green reading ofr october was incorrect.
% Tmatrix = [6.6635  -14.0425; -0.0267   32.1938];%% Calibrated 3/20/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated because new mobile screen was attached to the projector. Screen is now 1 inch farther than it used to be.
% Tmatrix = [6.3322  -13.4444; -0.0528   30.2447];%% Calibrated 3/22/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated because green projector was pushed up closer to the UV projector on the bretboard
Tmatrix = [6.3879  -13.5357; -0.0479   30.5975];%% Calibrated 3/26/18. Old EKB Projector, 405 LED. UV current 128, Green Current 128. 0.01" Teflon. Green 100% offset. Recalibrated because projectors were moved by allison, unclear what the movement was


if contrastNorm
    T=inv(Tmatrix); %Bring it back to the noninverted matrix
    alpha = sum(T(1,:))/sum(T(2,:)); %Ratio of the cone "responses" at the gray point.
    Tmatrix(:,2) = Tmatrix(:,2)/alpha; %Correct for the different gray points of the two cones
end

%Ugh, this is what I was doing prior to 5/14/16. It makes the S contrast
%too strong.  It also undersamples 90 to 180... the color directions
% Tmatrix(1,:) = Tmatrix(1,:)/norm(Tmatrix(1,:));
% Tmatrix(2,:) = Tmatrix(2,:)/norm(Tmatrix(2,:));

%This makes it so that the green and blue "guns" can't be > 1
normer = max([norm(Tmatrix(1,:)) norm(Tmatrix(2,:))]);
Tmatrix = Tmatrix/normer; %Norm of one row will be 1, and norm of the other row will be <1 

M = R*cos(theta*pi/180);  
S = R*sin(theta*pi/180);  
gb = Tmatrix*[M S]';
% greengain =  0.53*M + -.88*S;
% bluegain = 0*M + 2.9*S; 

