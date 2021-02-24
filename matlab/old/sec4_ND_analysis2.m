% solves the non-linear equation for lambda
clear;
% warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Set up directory (check if it exists)
foldername = 'sec4_ND_analysis5' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

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

% Specify delta:
delta = 0.030;
freq = solveOmegaND(delta, options);
Omega = freq.Omega;

% Obtain the eigenvalues for sampled dummy phases
phi1 = normrnd(0, delta, N, 1);
phi2 = normrnd(0, delta, N, 1);

Delta = bsxfun(@minus, phi1, phi2.');

% Solve eigenvalues
dist = solveEigsDetPlas(Omega, Delta, options);

% Information to export (all eigenvalues must scale up by scale*gain)
re = dist.re;
im = dist.im;
fvals = dist.fvals;
found = dist.found;
errors = dist.errors;
residual = dist.residual;
Delta = dist.Delta;
    
