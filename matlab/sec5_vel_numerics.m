% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec5_vel_numerics3' ;
subfoldername = 'plas'; % Change this to 'noplas'
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;
dir_subfolder = fullfile(cwd, 'data', foldername, subfoldername);

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder);
end

if ~exist(dir_subfolder, 'dir')
   mkdir(dir_subfolder);
end

% Parameters
N = 30;
omega0 = 1.0;
g = 1.5;
alphatau = 0.1; % Toggle this for plas, no plas
betatau = 0.005;
t0 = 0.0;
tf = 10000;
A = double(ones(N));
gain = 1.2; % 2.0;
vel0 = 1.0;
vel_min = 0.5;
vel_max = 150;
epsilon = 0.01;

asy = 0.05;

% Mean distance (increases)
distm_arr = 8.0; % 0.5:0.5:8.0;

parameters = struct('N', N, 'w0', omega0, 'g', g, 't0', t0, 'tf', tf, ...
    'vel0', vel0, 'vel_min', vel_min, 'vel_max', vel_max, ...
    'alphatau', alphatau, 'A', A, 'epsilon', epsilon, 'gain', gain, ...
    'betatau', betatau);

% History function (for each trial)
n_trials = numel(distm_arr);
std = 0.6; % Fix uniform bound [0,std] for init phases
L_freq = 0.25; % Multiple of g
freq_arr = rand(1,n_trials) - 0.5;
freq_arr = omega0 + L_freq * 2 * g * freq_arr;

% DDE options
options = ddeset();
options.NormControl = 'off';
% options.RelTol = 1e-6;
% options.AbsTol = 1e-6;
options.OutputFcn = @ddewbar;
options.MaxStep = 1.0;

% Wait bar
f = waitbar(0,'Starting trials...') ;

% MAIN LOOP
for k = 1:n_trials
    distm = distm_arr(k);
    freq = freq_arr(k);
        
    % Waitbar
    waittext = ['trial: ' num2str(k) ', distm = ' num2str(distm) ', freq = ' num2str(freq)] ;
    waitprog = k / n_trials ;
    waitbar(waitprog, f, waittext) ;

    % History function
    T = sqrt(3)*std;
    phases = T*rand(1,N) - T/2;
    Omega0 = freq;
    init_freqs = Omega0*ones(1,N);
    
    histfun = @(t) init_freqs*t + phases;
    % histob = IVPhistory(init_freqs, phases, parameters);
    % parameters.hist = histob.y;
    parameters.hist = histfun;
    
    % Adjust parameters
    dist = exprnd(distm, N);
    parameters.distm = distm;
    parameters.dist = dist;
    
    % Solve model
    sol = solveVelKuramotoPlas(parameters, options) ;
    
    % Export (transpose all matrices)
    t = sol.x.' ;
    y = sol.y(1:N,:).' ;
    yp = sol.yp(1:N,:).' ;

    vel = sol.y(N+1:end,:).' ;
    velp = sol.yp(N+1:end,:).' ;
    
    tau = (dist(:) ./ vel.').';
    
    hist_phases = phases.';
    hist_freq = Omega0;
    hist_std = std;
    
    % Asymptotic frequency and phases
    t_asy = t(t > (1-asy)*tf);
    y_asy = y(t > (1-asy)*tf,:);
    yp_asy = yp(t > (1-asy)*tf,:);
    
    t_diff = max(t_asy) - min(t_asy);
    
    Omega_asy = mean(trapz(t_asy, yp_asy)) / t_diff;
    
    % Save file
    filename = [num2str(k-1, '%02.f') '.mat'] ;
    dir_file = fullfile(dir_subfolder, filename) ;
    save(dir_file, 't', 'y', 'yp', 'vel', 'velp', 'tau', 'N', 'g', 'omega0', ...
        'dist', 'distm', 't0', 'tf', 'alphatau', 'hist_std', 'hist_freq', ...
        'hist_phases', 'Omega_asy', 'asy');
end

close(f);

% Save parameter file
filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_subfolder, filename0);
save(dir_file0, 'N', 'g', 'omega0', 'distm_arr', 't0', 'tf', 'alphatau', ...
    'vel0', 'vel_min', 'vel_max', 'gain', 'betatau', 'std', 'L_freq', ...
    'freq_arr', 'asy');

% TEST PLOT
% % plot(t, yp)

% figure;
% plot(t, sin(y - Omega_asy*t))
