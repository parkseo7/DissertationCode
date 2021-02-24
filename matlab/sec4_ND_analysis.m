% solves the non-linear equation for lambda
clear;
warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Set up directory (check if it exists)
foldername = 'sec4_ND_analysis4' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Parameters
gain = 80;
tau0 = 0.1;
g = 1.5;
omega0 = 1.0;

scale = 1.0; % Scale by z -> scale*gain*z

% Structure
options = struct('gain', gain, 'tau0', tau0, 'g', g, 'omega0', omega0);

% Obtain Omega(delta) array
deltaarr = 0.005:0.005:0.10; % 0.01:0.01:0.08; % Keep low
Omegaarr = omega0*ones(size(deltaarr));

for k = 1:numel(deltaarr)
    Omegaarr(k) = solveOmega(sqrt(2)*deltaarr(k), options);
end

% Eigenvalues at each synchronous point

% Numerical parameters
tol = 1.5; % Implement root finding on log10(err) < tol
roundingthreshold = 1;
imrange = [-100,100]; % Scaled range
rerange = [-80*g,80*g]; % Scaled range (actual range x scale*gain)
renum = 120;
imnum = 100;
Display = 'off';

MaxFunEvals = 300;
MaxIter = 300;

numberContours = 2^8;
MaxIntervalCount = 1e+6;

options = struct('tol', tol, 'rerange', rerange, 'imrange', imrange, ...
    'renum', renum, 'imnum', imnum, 'MaxFunEvals', MaxFunEvals, ...
    'MaxIter', MaxIter, 'roundingthreshold', roundingthreshold, ...
    'numberContours', numberContours, 'MaxIntervalCount', MaxIntervalCount, ...
    'scale', scale, 'Display', Display);

options.gain = gain;
options.tau0 = tau0;
options.g = g;
options.omega0 = omega0;

% For each sync point, obtain distribution and export
n_trials = numel(deltaarr);

% Wait bar
f = waitbar(0,'Starting trials...') ;
waitk = 0;

for i = 1:numel(deltaarr)
    
    % Sync state
    delta = deltaarr(i);
    Omega = Omegaarr(i);
    
    % Waitbar
    waittext = ['delta = ' num2str(delta), ' trial ' num2str(i) ' out of ' num2str(n_trials)];
    waitprog = i/ numel(deltaarr);
    waitbar(waitprog, f, waittext) ;
    
    % Solve eigenvalues
    dist = solveEigsND(Omega, sqrt(2)*delta, options);
    
    % Information to export (all eigenvalues must scale up by scale*gain)
    re = dist.re;
    im = dist.im;
    fvals = dist.fvals;
    found = dist.found;
    errors = dist.errors;
    residual = dist.residual;

    % Parameters and sync frequencies
    file_export = [num2str(i-1, '%02.f') '.mat']; 
    dir_file = fullfile(dir_folder, file_export);

    save(dir_file, 'g', 'gain', 'omega0', 'tau0', 'Omega', 'delta', ...
         're', 'im', 'fvals', 'found', 'errors', 'residual', ...
         'scale', 'tol')

end

close(f)

% Save parameter file
delta_arr = deltaarr.';
Omega_arr = Omegaarr.';
rerange = rerange.';
imrange = imrange.';
filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_folder, filename0);
save(dir_file0, 'g', 'omega0', 'gain', 'tau0', 'delta_arr', 'Omega_arr', ...
    'tol', 'scale', 'rerange', 'imrange', 'renum', 'imnum', 'MaxFunEvals', 'MaxIter');


