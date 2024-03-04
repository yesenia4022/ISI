function [theta, varargout] = circCorr(a, b, varargin)

%This comes from "Matteo Box"

% circCorr circular correlation coefficient based on Fisher & Lee Biometrika (1983).
%
% theta = circCorr(a, b) computes the circular correlation between vectors a
% and b. The units must be radians, and the data are assumed to be circular
% on a scale of 0-360 (i.e., the period is 2pi).
%
% theta = circCorr(a, b, period) let's you specify a period. Use pi to
% indicate that the data are circular on the scale of 0-180.
%
% [theta, p] = circCorr(a, b, period, p-ValueFlag) also returns a p-Value based
% on the Achieved Significance Level (ASL) as determined by permutation
% test.
%
% [theta, p, se] = circCorr(a, b, period, p-ValueFlag, seFlag) also returns the standard 
% error of theta as determined by bootstrapping 
%
% [theta, p, se, ci] = circCorr(a, b, period, p-ValueFlag, seFlag, ciFlag) also returns
% the 95% BCa bootstrap confidence interval
% 
% See also CIRCSTATS, CIRCSTATS360

% 2008-03 Steffen Katzner


if(nargin < 3)  %% set the period
    period = 2*pi;
    computePValue = 0;
    computeSe = 0;
    computeCi = 0;
elseif (nargin < 4)
    period = varargin{1};
    computePValue = 0;
    computeSe = 0;
    computeCi = 0;
elseif (nargin < 5)
    period = varargin{1};
    computePValue = varargin{2};
    computeSe = 0;
    computeCi = 0;
elseif (nargin < 6)
    period = varargin{1};
    computePValue = varargin{2};
    computeSe = varargin{3};
    computeCi = 0;
elseif (nargin < 7)
    period = varargin{1};
    computePValue = varargin{2};
    computeSe = varargin{3};
    computeCi = varargin{4};
end

%% Compute circular correlation coefficient
theta = cxcorr(a, b, period);

% figure;
% plot(a, b, 'bo');
% title(sprintf('theta = %2.3f', theta));


%% Compute p-value based on ASL

if computePValue
    v = [a, b]; % concatenate data
    g = [ones(1, length(a)), zeros(1, length(b))]; % indentity of data points
    B = 1000; % number of permutations

    gg = ones(B, length(g))*NaN; % each row represents one permutation
    thetaHat = ones(1, B)*NaN;
    for iPermutations = 1 : B
        rand('seed', iPermutations);
        [ignore, gg(iPermutations,:)] = sort(rand(1, length(g))); % random permutation of integers from 1 to length(g)
        g = g(gg(iPermutations,:)); % use random permutation to re-order g
        thetaHat(1, iPermutations) = cxcorr(v(g==1), v(g==0), period);
    end

    if size(gg, 1) ~= B
        error('permutations not unique');
    end

%     figure;
%     hist(thetaHat);

    % approximate ASL
    ind = find(thetaHat >= theta);
    p = length(ind) / B;
    
    varargout{1} = p;
end

%% Compute standard error

nBoot = 500; % number of bootstrap samples
if computeSe
    [bootstat, bootsam] = bootstrp(nBoot, @cxcorr, a, b, pi);
    se = std(bootstat);
    varargout{2} = se;
end

%% Compute 95% confidence interval

if computeCi
    ci = bootci(nBoot, @cxcorr, a, b, pi);
    varargout{3} = ci;
end






