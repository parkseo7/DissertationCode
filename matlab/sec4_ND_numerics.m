% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec4_ND_numerics1' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Implement stability analysis
is_stab = false;

% Parameters
N = 40;
omega0 = 1.0;
g = 1.5;
alphatau = 0.5;
t0 = 0.0;
tf = 500;
gain = 80;
tau0 = 0.1;
epsilon = 0.01;
A = double(ones(N));

parameters = struct('N', N, 'w0', omega0, 'g', g, 't0', t0, 'tf', tf, ...
    'gain', gain, 'tau0', tau0, 'epsilon', epsilon, 'alphatau', alphatau, 'A', A);

% Asymptotic values
asy = 0.1;

% Varying initial conditions
n_trials = 10; % Increase to 10

L_std = 0.5; % Uniform from [0, L_std]
L_freq = 0.25; % Multiple of g
std_arr = L_std*rand(1,n_trials);
freq_arr = rand(1,n_trials) - 0.5;
freq_arr = omega0 + L_freq * 2 * g * freq_arr;

% DDE options
options = ddeset();
options.NormControl = 'on';
% options.RelTol = 1e-6;
% options.AbsTol = 1e-6;
options.OutputFcn = @ddewbar;
options.MaxStep = 2.0 ;

% Eigenvalue search parameters (optional)
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

Eoptions = struct('logtol', logtol, 'rerange', rerange, 'imrange', imrange, ...
        'renum', renum, 'imnum', imnum, 'MaxFunEvals', MaxFunEvals, ...
        'MaxIter', MaxIter, 'roundingthreshold', roundingthreshold, ...
        'numberContours', numberContours, 'MaxIntervalCount', MaxIntervalCount, ...
        'cap', cap, 'scale', scale);

Eoptions.gain = gain;
Eoptions.tau0 = tau0;
Eoptions.g = g;
Eoptions.omega0 = omega0;
Eoptions.N = N;

% Wait bar
f = waitbar(0,'Starting trials...') ;

% MAIN LOOP
for k = 1:n_trials
    std = std_arr(k);
    freq = freq_arr(k);
        
    % Waitbar
    waittext = ['trial: ' num2str(k) ', std = ' num2str(std) ', freq = ' num2str(freq)] ;
    waitprog = k / n_trials ;
    waitbar(waitprog, f, waittext) ;

    % History function
    T = sqrt(3)*std;
    phases = T*rand(1,N) - T/2;
    Omega0 = freq;
    init_freqs = Omega0*ones(1,N);

    parameters.hist = IVPhistoryND(init_freqs, phases, tau0, parameters);

    % Solve model
    sol = solveDelayKuramotoPlas(parameters, options) ;
    
    % Export (transpose all matrices)
    t = sol.x.' ;
    y = sol.y(1:N,:).' ;
    yp = sol.yp(1:N,:).' ;

    tau = sol.y(N+1:end,:).' ;
    taup = sol.yp(N+1:end,:).' ;
    
    hist_phases = phases.';
    hist_freq = Omega0;
    hist_std = std;
    
    % Asymptotic frequency and phases
    t_asy = t(t > (1-asy)*tf);
    y_asy = y(t > (1-asy)*tf,:);
    yp_asy = yp(t > (1-asy)*tf,:);
    
    t_diff = max(t_asy) - min(t_asy);
    
    Omega_asy = mean(trapz(t_asy, yp_asy)) / t_diff;
    phases_asy = trapz(t_asy, y_asy - Omega_asy*t_asy) / t_diff;
    
    % Stability analysis
    eigs = 0;
    errors = 0;
    re = zeros(renum, imnum);
    im = zeros(renum, imnum);
    
    if is_stab
        Delta_asy = bsxfun(@minus, phases_asy, phases_asy.');
        dist = solveEigsDetPlas(Omega_asy, Delta_asy, Eoptions);
        eigs = dist.found;
        errors = dist.errors;
        re = dist.re;
        im = dist.im;
    end
    
    % Save file
    filename = [num2str(k-1, '%02.f') '.mat'] ;
    dir_file = fullfile(dir_folder, filename) ;
    save(dir_file, 't', 'y', 'yp', 'tau', 'taup', 'N', 'g', 'omega0', ...
        'tau0', 'gain', 't0', 'tf', 'alphatau', 'hist_std', 'hist_freq', 'hist_phases', ...
        'Omega_asy', 'phases_asy', 'asy', 'eigs', 'errors', 're', 'im');
end

close(f);

% Save parameter file
rerange = rerange.';
imrange = imrange.';
filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_folder, filename0);
save(dir_file0, 'N', 'g', 'omega0', 'tau0', 'gain', 't0', 'tf', 'alphatau', ...
    'epsilon', 'logtol', 'rerange', 'imrange', 'cap', 'renum', 'imnum', ...
    'MaxFunEvals', 'MaxIter', 'L_std', 'L_freq', 'asy');