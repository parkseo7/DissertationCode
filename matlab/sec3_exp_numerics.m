% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec3_exp_numerics4' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% DETERMINE IF SEARCH FOR THEORETICAL EIGENVALUES
is_eigs = true;

% FIXED PARAMETERS
N = 60;
g = 1.5;
omega0 = 1.0;
t0 = 0;
tf = 500;
A = ones(N);

% Structure
parameters = struct('g', g, 'w0', omega0, 'tau0', 0, 'A', A, 't0', t0, 'tf', tf);

% History function conditions
std = 2.0;
T = sqrt(3)*std;
phases = T*rand(1,N) - T/2;
Omega0 = omega0;
init_freqs = Omega0*ones(1,N);
phi0 = phases;

% VARIABLE PARAMETERS
taum_arr = 0:0.5:7 ; % linspace(0, 7, 21);

eqOmega = omega0;

% EIGENVALUE OPTIONS

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

found = 0;
errors = 0;

% DDE options
ddeopts = ddeset() ;
% ddeopts.NormControl = 'on';
ddeopts.OutputFcn = @ddewbar;
ddeopts.MaxStep = 2.0 ;

% Wait bar
f = waitbar(0,'Starting trials...') ;
n_trials = numel(taum_arr);

% MAIN LOOP
for j = 1:n_trials
    taum = taum_arr(j);

    % Waitbar
    waittext = ['taum = ' num2str(taum), ', num ' num2str(j) ' out of ' num2str(n_trials)] ;
    waitprog = j / n_trials ;
    waitbar(waitprog, f, waittext) ;

    % History function
    parameters.tau0 = taum; % Use the mean delay to establish history
    parameters.hist = IVPhistory(init_freqs, phases, parameters);

    % Delays
    parameters.g = g;
    parameters.tau = exprnd(taum, N, N);
    tau = parameters.tau;
    
    % Theoretical sync. frequency
    if taum > 0
        mu = 1 / taum;
    else
        mu = 1;
    end
    p = [1 -omega0 mu^2+g*mu -omega0*mu^2];
    eqOmegas = roots(p).';
    
    % Keep the real root (if any)
    for k=1:numel(eqOmegas)
        if isreal(eqOmegas(k))
            eqOmega = real(eqOmegas(k));
        end
    end
    
    % Solve eigenvalues (cubic)
    q1 = 1;
    q2 = 2*mu + g*mu^2/(mu^2+eqOmega^2);
    q3 = (mu^2 + eqOmega^2) + 2*g*mu^3/(mu^2 + eqOmega^2) - g*mu;
    q4 = 0;
    q = [q1 q2 q3 q4];
    eigenvalues = roots(q).';
    
    % Eigenvalue distribution
    options.tau0 = parameters.tau;
    options.taum = taum;
    
    % Get distribution
    if is_eigs
        dist = solveEigsDet(eqOmega, options);

        % Set up values
        found = dist.found;
        errors = dist.errors;
    end
    
    % Solve model
    sol = solveDelayKuramoto(parameters, ddeopts) ;
    
    % Export (transpose all matrices)
    t = sol.x.' ;
    y = sol.y(1:N,:).' ;
    yp = sol.yp(1:N,:).' ;

    % Save file
    filename = [num2str(j-1, '%02.f') '.mat'] ;
    dir_file = fullfile(dir_folder, filename) ;
    save(dir_file, 't', 'y', 'yp', 'N', 'omega0', 'g', 'taum', ...
        'tau', 'phi0', 't0', 'tf', 'std', 'Omega0', 'eqOmegas', ...
        'eqOmega', 'eigenvalues', 'found', 'errors')
end

close(f);

% Save parameter file
taum = taum_arr.';

filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_folder, filename0);
save(dir_file0, 'N', 'omega0', 'g', 'taum', 'phi0', 't0', 'tf', 'std', 'Omega0');