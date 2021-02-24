% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec4_ND_injury' ;
subfoldername = 'noplas'; % Change this to 'noplas'
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
N = 40;
omega0 = 1.0;
g = 1.5;
alphatau = 0.5;
t0 = 0.0;
tf = 500;
gain = 0; % Toggle this to 0 with no plasticity
tau0 = 2.0;
epsilon = 0.01;

inj = 0.0;
tinjs = 240;
tinjf = 260;
    
parameters = struct('N', N, 'w0', omega0, 'g', g, 't0', t0, 'tf', tf, ...
    'gain', gain, 'tau0', tau0, 'epsilon', epsilon, 'alphatau', alphatau, ...
    'inj', inj, 'tinjs', tinjs, 'tinjf', tinjf, 'A', ones(N));

% Asymptotic values
asy = 0.1;

% Fixed initial condition
std = 0.5;
initfreq = 0.8;

% Vary injury level
inj_arr = 0.0:0.1:0.9; % 0.0:0.1:1.0
n_trials = numel(inj_arr);

% DDE options
options = ddeset();
options.NormControl = 'on'; % Keep on
options.OutpuFcn = @ddewbar;
% options.MaxStep = 1.0 ;

% Wait bar
f = waitbar(0,'Starting trials...') ;

% MAIN LOOP
for k = 1:n_trials
    
    inj = inj_arr(k);
    parameters.inj = inj;
    
    % Waitbar
    waittext = ['trial: ' num2str(k) ', inj = ' num2str(inj)] ;
    waitprog = k / n_trials ;
    waitbar(waitprog, f, waittext) ;

    % History function
    T = sqrt(3)*std;
    phases = T*rand(1,N) - T/2;
    initfreqs = initfreq*ones(1,N);

    parameters.hist = IVPhistory(initfreqs, phases, parameters);

    % Solve model
    sol = solveDelayKuramotoInj(parameters, options) ;
    
    % Export (transpose all matrices)
    t = sol.x.' ;
    y = sol.y(1:N,:).' ;
    yp = sol.yp(1:N,:).' ;
    
    tau = sol.y(N+1:end,:).' ;
    taup = sol.yp(N+1:end,:).' ;
    
    hist_phases = phases.';
    
    A_inj = sol.A_inj;
    molarr = sol.molarr;
    
    % Asymptotic frequency and phases
    tinj_post = tf - (tf - tinjf)*asy;
    tinj_pre = tinjs - (tinjs - t0)*asy;
   
    tasy_pre = t((tinj_pre < t) & (t < tinjs));
    yasy_pre = y((tinj_pre < t) & (t < tinjs),:);
    ypasy_pre = yp((tinj_pre < t) & (t < tinjs),:);
    tdiff_pre = max(tasy_pre) - min(tasy_pre);
    
    tasy_post = t(t > tinj_post);
    yasy_post = y(t > tinj_post,:);
    ypasy_post = yp(t > tinj_post,:);
    tdiff_post = max(tasy_post) - min(tasy_post);
    
    % MUST TAKE ASYMPTOTIC PRE-INJURY, POST-INJURY
    Omega_pre = mean(trapz(tasy_pre, ypasy_pre)) / tdiff_pre;
    phases_pre = trapz(tasy_pre, yasy_pre - Omega_pre*tasy_pre) / tdiff_pre;
    
    Omega_post = mean(trapz(tasy_post, ypasy_post)) / tdiff_post;
    phases_post = trapz(tasy_post, yasy_post - Omega_post*tasy_post) / tdiff_post;
    
    % Save file
    filename = [num2str(k-1, '%02.f') '.mat'] ;
    dir_file = fullfile(dir_subfolder, filename) ;
    save(dir_file, 't', 'y', 'yp', 'tau', 'taup', 'N', 'g', 'omega0', ...
        'tau0', 'gain', 't0', 'tf', 'alphatau', 'inj', 'tinjs', 'tinjf', ...
        'hist_phases', 'A_inj', 'molarr', 'Omega_pre', 'phases_pre', ...
        'Omega_post', 'phases_post', 'asy');
end

close(f);

hist_freq = initfreq;
hist_std = std;
    
% Save parameter file
filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_subfolder, filename0);
save(dir_file0, 'N', 'g', 'omega0', 'tau0', 'gain', 't0', 'tf', 'alphatau', ...
    'epsilon', 'hist_std', 'hist_freq', 'inj_arr', 'tinjs', 'tinjf', 'asy');