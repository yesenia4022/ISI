function s_mu = ktracker(z,T)

% s_mu =  [z(:,1); (z(:,2)-z(:,1))/T; 2*((z(:,3)-z(:,2))-(z(:,2)-z(:,1)))/(T^2)];
% D = shat*shat';
% P = inv(inv(D) + H'*inv(R)*H);

Ord = 3;

N = length(z(1,:));     %Samples
D = length(z(:,1));     %Dimensions
V = Ord*D;              %Length of State Vector

Qn = .1;
Q = eye(V);  
Q(1:D,:) = Qn*Q(1:D,:);
Q(D+1:2*D,:) = .02*Qn*Q(D+1:2*D,:);
Q(2*D+1:3*D,:) = .005*Qn*Q(2*D+1:3*D,:);  %covariance of the process noise (lower values smooth more)

Rn = 100;  
R = Rn*eye(D);  %covariance of the observation noise (higher values smooth more)
Ri = inv(R);
H = eye(V);
H = H(:,1:D)';

F = [1 zeros(1,D-1) T zeros(1,D-1) .5*T^2 zeros(1,D-1)];

for i = 1:V-1
    Fnext = [zeros(1,i) F(1,1:V-i)];
    F = [F;Fnext];
end

%INITIAL CONDITIONS:
z_hat = z(:,1);
P_mu = 0*eye(V);
s_mu = [z(:,1); zeros(V-D,1)];
e = z_hat(:,1);

%KALMAN FILTER RECURSION (REAL-TIME LOOP):

for cycle = 1:N-1
    
    P_tu = F*P_mu*F' + Q; 
    Re = R + H*P_tu*H';
    P_mu = P_tu - P_tu*H'*inv(Re)*H*P_tu;
    
    s_tu = F*s_mu(:,cycle);
    z_hat(:,cycle+1) = H*s_tu;   %time_update of state

    e(:,cycle+1) = z(:,cycle+1) - H*s_tu;
    s_mu(:,cycle+1) = s_tu + P_mu*H'*Ri*e(:,cycle+1);
    
end

%DETERMINE CONFIDENCE IN THE WHITENESS:
% N = length(e(1,:));
% Cx = (1/N)*conv(fliplr(e(1,:)),e(1,:));
% Cx = Cx(N:length(Cx));
% rhox = Cx/Cx(1);
% Confx = (sum(sign(abs(rhox(2:N)) - (1.96/(N^.5))) + 1))/2;
% Confx = (1 - Confx/N)*100;
% Cy = (1/N)*conv(fliplr(e(2,:)),e(2,:));
% Cy = Cy(N:length(Cy));
% rhoy = Cy/Cy(1);
% Confy = (sum(sign(abs(rhoy(2:N)) - (1.96/(N^.5))) + 1))/2;
% Confy = (1 - Confy/N)*100;


% x_hat = z_hat(1,:);
% y_hat = z_hat(2,:);
% x_error = x - x_hat;
% y_error = y - y_hat;


