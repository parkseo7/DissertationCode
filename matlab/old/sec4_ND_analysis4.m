% solves the non-linear equation for lambda
clear;
% warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Parameters
N = 40; % 40;
gain = 80;
tau0 = 0.1;
g = 1.5;
omega0 = 1.0;

% Numerical parameters
logtol = -5; % Implement root finding on all log10err < tol
scale = gain; % Scales lambda by gain*lambda
rerange = [-g,g]; % Scaled re range

imrange = [-50,50]; % Scaled im range
cap = 1; % Cap for det (found only if < cap)
renum = 60;
imnum = 100;
MaxFunEvals = 2000;
MaxIter = 2000;

roundingthreshold = 6;
numberContours = 2^8;
MaxIntervalCount = 1e+6;

options = struct('logtol', logtol, 'rerange', rerange, 'imrange', imrange, ...
    'renum', renum, 'imnum', imnum, 'MaxFunEvals', MaxFunEvals, ...
    'MaxIter', MaxIter, 'roundingthreshold', roundingthreshold, ...
    'numberContours', numberContours, 'MaxIntervalCount', MaxIntervalCount, ...
    'cap', cap, 'scale', scale);

% Add parameters
options.gain = gain;
options.tau0 = tau0;
options.g = g;
options.omega0 = omega0;
options.N = N;
options.get_arr = false;

% Iteration settings
options.num_iter = 1000;
options.learning = 0.1;

% Specify delta:
delta = 0.050;
freq = solveOmegaND(delta, options);
Omega = freq.Omega;

% Initial phases
phi = normrnd(0, delta, 1, N);

% Optimize phases
dist = solvePhis(Omega, phi, options);
phi_c = dist.phi;

% TEST VARIABLES
Delta = bsxfun(@minus,phi_c,phi_c');
sinDelta = sin(-Omega*max(tau0 + gain*sin(Delta), 0) + Delta);
SUM1 = sum(sinDelta, 1);
SUM2 = sum(sinDelta, 2);

% disp(Omega - omega0 - (g/N)*SUM2);
disp(sum((Omega - omega0 - (g/N)*SUM2).^2));