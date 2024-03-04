function [S_rod,S_Mcone,S_Scone, lambda] = photoreceptor_templates(rod_peak,M_peak,S_peak) 

% rod_peak ~ 500;
% M_peak ~ 508;
% S_peak ~ 360;

%From Govardovskii et al 2000...

lambda = 350:700;
A = 69.7;
B = 28;
b = .922;
C = -14.9;
c = 1.104;
D = 0.674;
A_beta = 0.26;

%% Rod
lambda_peak = rod_peak;
x = 1./(lambda/lambda_peak);
a = .8795 + .0459*exp(-(lambda_peak-300)^2/11940);
lambda_beta = 189 + .315*lambda_peak;
d = -40.5+.195*lambda_peak;

S_alpha = [exp(A*(a-x)) + exp(B*(b-x)) + exp(C*(c-x)) + D].^-1;
S_beta = A_beta*exp( -((lambda-lambda_beta)/d).^2);
S_rod = S_alpha + S_beta;

% figure, plot(lambda,S_rod,'k')

%% M cone
lambda_peak = M_peak;
x = 1./(lambda/lambda_peak);
a = .8795 + .0459*exp(-(lambda_peak-300)^2/11940);
lambda_beta = 189 + .315*lambda_peak;
d = -40.5+.195*lambda_peak;

S_alpha = [exp(A*(a-x)) + exp(B*(b-x)) + exp(C*(c-x)) + D].^-1;
S_beta = A_beta*exp( -((lambda-lambda_beta)/d).^2);
S_Mcone = S_alpha + S_beta;

% hold on, plot(lambda,S_Mcone,'g')
% 
% legend('rod','M cone')

%% S cone
lambda_peak = S_peak;
x = 1./(lambda/lambda_peak);
a = .8795 + .0459*exp(-(lambda_peak-300)^2/11940);
lambda_beta = 189 + .315*lambda_peak;
d = -40.5+.195*lambda_peak;

S_alpha = [exp(A*(a-x)) + exp(B*(b-x)) + exp(C*(c-x)) + D].^-1;
S_beta = A_beta*exp( -((lambda-lambda_beta)/d).^2);
S_Scone = S_alpha + S_beta;

% hold on, plot(lambda,S_Scone,'b')
% 
% legend('rod','M cone','S cone')
