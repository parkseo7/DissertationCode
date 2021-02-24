clear;
warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Set up directory (check if it exists)
foldername = 'sec3_exp_analysis' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Fixed parameters
N = 60;
g = 1.5;
omega0 = 1.0;

% Numerical parameters
tol = -10; % Implement root finding on all log10err < tol
rerange = [-g,g];
imrange = [-100,100];
cap = 1; % Cap for det (found only if < cap)
renum = 60;
imnum = 100;
MaxFunEvals = 2000;
MaxIter = 2000;

roundingthreshold = 6;
numberContours = 2^8;
MaxIntervalCount = 1e+6;

options = struct('tol', tol, 'rerange', rerange, 'imrange', imrange, ...
    'renum', renum, 'imnum', imnum, 'MaxFunEvals', MaxFunEvals, ...
    'MaxIter', MaxIter, 'roundingthreshold', roundingthreshold, ...
    'numberContours', numberContours, 'MaxIntervalCount', MaxIntervalCount, ...
    'cap', cap);

options.g = g;
options.omega0 = omega0;
options.N = N;

% Variable parameters
taum_step = 0.5;
taum_end = 8;
taum_arr = taum_step:taum_step:taum_end;

% MAIN LOOP:

% Wait bar
f = waitbar(0,'Starting trials...') ;

num_trials = numel(taum_arr);

% MAIN LOOP
for j = 1:n_trials
    taum = taum_arr(j);

    % Waitbar
    waittext = ['taum = ' num2str(taum), ', num ' num2str(j) ' out of ' num2str(n_trials)] ;
    waitprog = j / n_trials ;
    waitbar(waitprog, f, waittext) ;
    
    % Theoretical sync. frequency
    mu = 1 / taum;
    p = [1 -omega0 mu^2+g*mu -omega0*mu^2];
    eqOmegas = roots(p).';
    
    % Keep the real root (if any)
    for k=1:numel(eqOmegas)
        if isreal(eqOmegas(k))
            eqOmega = real(eqOmegas(k));
        end
    end
    
    % Solve expected eigenvalues
    q1 = 1;
    q2 = 2*mu + g*mu^2/(mu^2+eqOmega^2);
    q3 = (mu^2 + eqOmega^2) + 2*g*mu^3/(mu^2 + eqOmega^2) - g*mu;
    q4 = 0;
    q = [q1 q2 q3 q4];
    eigenvalues = roots(q);
    
    tau0 = exprnd(taum, N);
    options.tau0 = tau0;
    options.taum = taum;
    
    % Get distribution
    dist = solveEigsDet(eqOmega, options);
    
    % Set up values
    found = dist.found;
    errors = dist.errors;
    
    re = dist.re;
    im = dist.im;
    fvals = dist.fvals;
    zero = dist.zero;
    
    % Save file
    filename = [num2str(j-1, '%02.f') '.mat'] ;
    dir_file = fullfile(dir_folder, filename) ;
    save(dir_file, 'N', 'g', 'omega0', 'taum', 'tau0', 'eqOmega', ...
        'eigenvalues', 'found', 'errors', 're', 'im', 'fvals', 'zero');
end

close(f);

% Save parameter file
taum = taum_arr.';
rerange = rerange.';
imrange = imrange.';
filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_folder, filename0);
save(dir_file0, 'N', 'g', 'omega0', 'taum', 'tol', 'rerange', 'imrange', ...
    'cap', 'renum', 'imnum', 'MaxFunEvals', 'MaxIter');

% figure();
% hold on;
% scatter(real(dist.found), imag(dist.found))
% scatter(real(eigenvalues), imag(eigenvalues))
% title(['$$N = ' num2str(N) ', \Omega = ' num2str(eqOmega) ', \tau_m = ' num2str(taum) '$$'],'interpreter','latex')
    

