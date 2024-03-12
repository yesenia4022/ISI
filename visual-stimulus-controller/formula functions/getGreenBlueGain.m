function gb = getGreenBlueGain(theta,R,contrastNorm)

global Mstate

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


%Tmatrix =  [16.0138  -7.1109; -0.0273  10.2242]; %Calibrated 6/24/17. Moved projectors closer .01" Teflon (30" from screen). New UV projector (EKB) with 385 LED. 0% offset RGB.  Green at 64, UV at 100

% Tmatrix =    [16.6250  -40.2109; -0.0509   35.9227]; %Calibrated 7/13/17. .01" Teflon. Proj 30" from screen. Newer UV projector (EKB) with 410 LED. 0% offset RGB.  Green at 64, UV at 128 current

% Tmatrix = [4.9660  -11.8357; -0.0181   19.4845]; % Calibrated 9/20/17. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration. Green at 128, UV at 128 current
% Tmatrix = [5.0187  -12.0015; -0.0269   20.0749]; % Calibrated 1/2/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration. Green at 128, UV at 128 current
% Tmatrix = [1.7624011695577         -4.52542124192665; -0.00312975128997142            7.189287216622]; % % Calibrated 1/5/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration. Green at 244, UV at 250 current

% Tmatrix = [3.43823698791844 -9.01330508560316; -0.0207073654009068 15.0763254850175]; % Calibrated 1/15/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration. Green at 128, UV at 128 current
% Tmatrix = [1.76878808928921 -4.65714016547708; -0.00327750006169761 7.63286855936557]; % Calibrated 1/15/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration. Green at 244, UV at 250 current

% Tmatrix = [5.1322 -11.9281; -0.0221 22.2024]; % Calibrated 2/26/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration but UV LED swapped for 'UV405-keynote' from the Keynote (Dallas?) projector. Green at 128, UV at 128 current
% Tmatrix = [2.9201 -6.7603; -0.0091 12.8109]; % Calibrated 2/26/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration but UV LED swapped for 'UV405-keynote' from the Keynote (Dallas?) projector. Green at 220, UV at 220 current
% Tmatrix = [2.6785 -6.2333; -0.0097 11.4866]; % Calibrated 2/26/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration but UV LED swapped for 'UV405-keynote' from the Keynote (Dallas?) projector. Green at 240, UV at 245 current
% Tmatrix = [2.9201 -6.7878; -0.0089 12.5268]; % Calibrated 2/26/18. .01" Teflon, 43cm from screen, same UV and RGB projectors from 7/13/17 calibration but UV LED swapped for 'UV405-keynote' from the Keynote (Dallas?) projector. Green at 220, UV at 225 current

% Tmatrix = [8.4109  -18.9026; -0.0428   25.0859]; % Calibrated 4/3/18. .01" Teflon 54cm from screen, new nUV LED (EKB) projector, same RGB projector. (G128, UV128)

% Tmatrix = [13.1343  -28.4886; -0.1583  174.0824]; % Calibrated 4/10/2018.Same projector setup as directly above, but using third of exposure time. (G128, UV128)

%Tmatrix = [4.6506  -10.4564; -0.0087   13.9312]; % Calibrated 4/3/18. .01" Teflon 54cm from screen, new nUV LED (EKB) projector, same RGB projector. (G223, UV225)


Tinv = inv(Mstate.T_MSGU(1:2,:));  %inverse of outer product of cones x gun spectra: [Mcone Scone]'*[Green UV];)

if contrastNorm
    T=inv(Tinv); %Bring it back to the noninverted matrix
    alpha = sum(T(1,:))/sum(T(2,:)); %Ratio of the cone "responses" at the gray point.
    Tinv(:,2) = Tinv(:,2)/alpha; %Correct for the different gray points of the two cones
end

%Ugh, this is what I was doing prior to 5/14/16. It makes the S contrast
%too strong.  It also undersamples 90 to 180... the color directions
% Tmatrix(1,:) = Tmatrix(1,:)/norm(Tmatrix(1,:));
% Tmatrix(2,:) = Tmatrix(2,:)/norm(Tmatrix(2,:));

%This makes it so that the green and blue "guns" can't be > 1
normer = max([norm(Tinv(1,:)) norm(Tinv(2,:))]);
Tinv = Tinv/normer; %Norm of one row will be 1, and norm of the other row will be <1 

M = R*cos(theta*pi/180);
S = R*sin(theta*pi/180);
gb = Tinv*[M S]';
% greengain =  0.53*M + -.88*S;
% bluegain = 0*M + 2.9*S; 

