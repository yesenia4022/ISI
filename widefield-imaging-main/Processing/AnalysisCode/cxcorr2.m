function r = cxcorr2(a,b,varargin)

%2 uses a "coherence" (vector strength) type of metric to gauge whether the pairs are more or
%less similar
%This is quite different than the regular cxcorr... its based on the phase
%shift instead of the slope of the scatter plot.

if(nargin>2)  %% set the period
    per = varargin{1};
else
    per = 2*pi;
end

a = a/per*2*pi;
b = b/per*2*pi;

coh = mean(exp(1i*a).*exp(-1i*b)); %coherency
r = abs(coh)*cos(angle(coh));  %projection of coherency vector

%%
%Dario's doesn't go to zero at a 90 deg phase shift... seems biased
%positive in general
%r = mean(abs(exp(1i*a) + exp(1i*b)) - 1); %Paik and Ringach 2011

