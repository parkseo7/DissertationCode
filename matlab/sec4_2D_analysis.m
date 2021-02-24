% Script to export function and error meshes of 2D Kuramoto analysis
% to be plotted.
clear;
warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Set up directory (check if it exists)
foldername = 'sec4_2D_analysis1' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Parameters
gain = 30;
tau0 = 0.5;
g = 1.5/2;
omega0 = 1.0;

parameters = struct('gain', gain, 'tau0', tau0, 'g', g, 'omega0', omega0);

% Options for finding synchronization states
options1 = struct();
options1.numsteps = 2000;
options1.roundingthreshold = 4;

% Synchronization states
sync = solveOmega2D(parameters, options1);

% Export parameters and sync frequencies
file_sync = 'sync_freqs.mat';
dir_file0 = fullfile(dir_folder, file_sync);

Omegas = sync.Omega;
Deltas = sync.Delta;
Omega_arr = sync.x_arr;
error_arr = sync.y_arr;
error = sync.error;

save(dir_file0, 'g', 'gain', 'omega0', 'tau0', 'Omegas', 'Deltas', 'Omega_arr',...
    'error_arr', 'error');

% Options for finding eigenvalues
options2 = struct();
options2.tol = 1e0;
options2.roundingthreshold = 3;

options2.rerange = 3;
options2.imrange = 6;

options2.renum = 50;
options2.imnum = 50;

% fminsearch options to find zeroes
fminOptions = struct;
fminOptions.MaxIter = 2000;
fminOptions.MaxFunEvals = 1000;

% Eigenvalues (for all sync frequencies < w0)
Omegainds = find(Omegas < omega0);
Omegas = Omegas(Omegainds);
Deltas = Deltas(Omegainds);

% Wait bar
f = waitbar(0,'Starting trials...') ;
n_trials = numel(Omegas);

for i = 1:n_trials
    Omega = Omegas(i);
    Delta = Deltas(i);
    
    dist = solveEigs2D(Omega, Delta, parameters, options2, fminOptions);
    
    % Waitbar
    waittext = ['sync Omega = ' num2str(Omega), ', num ' num2str(i) ' out of ' num2str(n_trials)] ;
    waitprog = i / n_trials;
    waitbar(waitprog, f, waittext);
    
    % Set parameters
    eigs = dist.found;
    eigerrors = dist.foundVals;
    polyroots = dist.polyroots;
    re = dist.re;
    im = dist.im;
    fvals = dist.fvals;
    
    % Save
    file_eigs = ['eig_dist' num2str(i) '.mat'];
    dir_file1 = fullfile(dir_folder, file_eigs);
    save(dir_file1, 'Omega', 'Delta', 'eigs', 'eigerrors', 'polyroots', 're', 'im', 'fvals')

end

close(f)



