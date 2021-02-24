% Set up directory (check if it exists)
foldername = 'sec4_2D_4' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Add to function path
addpath('fcns');

% Constant parameters
par = struct ;
par.w0 = 1.0 ;
par.g = 1.5 ;
par.t0 = 0 ;
par.tf = 15000 ;
par.gain = 30 ;
par.alphatau = 0.001 ;
par.tau0 = 0.1;
par.epsilon = 0.01;
par.A = [0 1; 1 0];
par.offset = 1.0;

% Export names
g = par.g ;
gain = par.gain ;
omega0 = par.w0 ;
tau0 = par.tau0.' ;
alphatau = par.alphatau;
t0 = par.t0 ;
tf = par.tf ;
offset = par.offset;

% Varying parameters
n_Delta = 4;
n_freq = 5;
L_Delta = 1.0;
L_freq = 0.5;
Delta_arr = L_Delta * linspace(0,1, n_Delta+1);
Delta_arr = Delta_arr(2:end);

freq_arr = linspace(0,1, n_freq) - 0.5;
freq_arr = par.w0 + L_freq * 2 * g * freq_arr;

total = n_Delta * n_freq;

% DDE options
ddeopts = ddeset() ;
% ddeopts.MaxStep = 1.0 ;

% Wait bar
f = waitbar(0,'Starting trials...') ;

% MAIN LOOP
for k = 1:n_Delta
    for l = 1:n_freq
       
        Delta0 = Delta_arr(k);
        init_freq = freq_arr(l);

        % Waitbar
        trial_num = (k-1)*n_freq + l;
        waittext = ['Delta = ' num2str(Delta0) ', Init.freq = ' num2str(init_freq)] ;
        waitprog = trial_num / n_trials;
        waitbar(waitprog, f, waittext) ;

        % Solve model
        par.hist = IVPhistory2D(init_freq, Delta0, par);
        % par.hist = @(t) t * [init_freq, init_freq] + [0, Delta0];
        
        sol = solvemodel2D(par, ddeopts) ;

        % Export (transpose all matrices)
        t = sol.x.' ;
        y = sol.y(1:2,:).' ;
        yp = sol.yp(1:2,:).' ;

        tau = sol.y(3:end,:).' ;
        taup = sol.yp(3:end,:).' ;

        % Save file
        filename = ['2D_num_' num2str(trial_num) '.mat'] ;
        dir_file = fullfile(dir_folder, filename) ;
        save(dir_file, 't', 'y', 'yp', 'tau', 'taup', 'tau0', 'gain', 'omega0', ...
            'g', 'tf', 'Delta0', 'init_freq', 'alphatau', 'offset')
    end
end

close(f)