% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec3_dirac_numerics6' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% FIXED PARAMETERS
N = 20;
omega0 = 1.0;
t0 = 0;
tf = 200;
A = ones(N);

% Structure
parameters = struct('g', 0, 'w0', omega0, 'tau0', 0, 'A', A, 't0', t0, 'tf', tf);

% History function conditions
std = pi/2;
T = sqrt(3)*std;
phases = T*rand(1,N) - T/2;
Omega0 = omega0;
init_freqs = Omega0*ones(1,N);
phi0 = phases;

% VARIABLE PARAMETERS
tau0_arr = 0:0.25:12;
g_arr = 0:0.1:1.5;

% DDE options
ddeopts = ddeset() ;
% ddeopts.NormControl = 'on';
ddeopts.OutputFcn = @ddewbar;
% ddeopts.MaxStep = 3.0 ;

% Wait bar
f = waitbar(0,'Starting trials...') ;

n_trials = numel(tau0_arr)*numel(g_arr);
n_g = numel(g_arr);

% MAIN LOOP
for i = 1:numel(g_arr)
    for j = 1:numel(tau0_arr)
        tau0 = tau0_arr(j);
        g = g_arr(i);
            
        % Waitbar
        trial_num = (i-1)*n_g + j;
        waittext = ['g = ' num2str(g), ', tau0 = ' num2str(tau0)] ;
        waitprog = trial_num / n_trials ;
        waitbar(waitprog, f, waittext) ;

        % History function
        parameters.tau0 = tau0;
        parameters.hist = IVPhistory(init_freqs, phases, parameters);
    
        % Delays
        parameters.g = g;
        parameters.tau = tau0*ones(N);

        % Solve model
        sol = solveDelayKuramoto(parameters, ddeopts) ;

        % Export (transpose all matrices)
        t = sol.x.' ;
        y = sol.y(1:N,:).' ;
        yp = sol.yp(1:N,:).' ;

        % Save file
        filename = [num2str(i-1, '%02.f') '_' num2str(j-1, '%02.f') '.mat'] ;
        dir_file = fullfile(dir_folder, filename) ;
        save(dir_file, 't', 'y', 'yp', 'N', 'omega0', 'g', ...
            'tau0', 'phi0', 't0', 'tf', 'std', 'Omega0')
    end
end

close(f);

% Save parameter file
g = g_arr.';
tau0 = tau0_arr.';

filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_folder, filename0);
save(dir_file0, 'N', 'omega0', 'g', 'tau0', 'phi0', 't0', 'tf', 'std', 'Omega0');

