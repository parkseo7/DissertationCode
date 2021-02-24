% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec4_2D_numerics' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% FIXED PARAMETERS
par = struct ;
par.t0 = 0 ;
par.w0 = 1.0 ;
par.g = 1.5/2 ;
par.gain = 30 ;
par.tau0 = 2.0;
par.alphatau = 0.5 ; % Toggle this
par.tf = round(50*par.alphatau^(-1), 0) ; % Toggle this
par.epsilon = 0.01; 
par.offset = 0.0;

% Export names
g = par.g ;
gain = par.gain ;
omega0 = par.w0 ;
tau0 = par.tau0.' ;
alphatau = par.alphatau;
t0 = par.t0 ;
tf = par.tf ;
offset = par.offset;

% Varying initial conditions
n_Omega = 7;
n_Delta = 11;
Omega0_arr = linspace(omega0-g, omega0+g, n_Omega);
Delta0_arr = linspace(0, 1, n_Delta);

n_trials = n_Omega * n_Delta;

% DDE options
options = ddeset() ;
options.NormControl = 'on';
options.MaxStep = 1.0 ;

% Wait bar
f = waitbar(0,'Starting trials...') ;


% MAIN LOOP
for i = 1:n_Delta
    for j = 1:n_Omega
        Omega0 = Omega0_arr(j);
        Delta0 = Delta0_arr(i);
            
        % Waitbar
        trial_num = (i-1)*n_Delta + j;
        waittext = ['Omega0 = ' num2str(Omega0), ', Delta0 = ' num2str(Delta0)] ;
        waitprog = trial_num / n_trials ;
        waitbar(waitprog, f, waittext) ;
        
        % Solve model
        par.hist = IVPhistory2D(Omega0, Delta0, par);
        sol = solvemodel2D(par, options) ;

        % Export (transpose all matrices)
        t = sol.x.' ;
        y = sol.y(1:2,:).' ;
        yp = sol.yp(1:2,:).' ;

        tau = sol.y(3:end,:).' ;
        taup = sol.yp(3:end,:).' ;
        
        % Save file
        filename = [num2str(i-1, '%02.f') '_' num2str(j-1, '%02.f') '.mat'] ;
        dir_file = fullfile(dir_folder, filename) ;
        save(dir_file, 't', 'y', 'yp', 'tau', 'taup', 'Omega0', 'Delta0', ... 
             'tau0', 'gain', 'omega0', 'g', 't0', 'tf', 'alphatau', 'offset')
    end
    
end

close(f);

% Save parameter file
Omega0 = Omega0_arr.';
Delta0 = Delta0_arr.';

filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_folder, filename0);
save(dir_file0, 'tau0', 'gain', 'omega0', 'g', 't0', 'tf', ... 
     'Delta0', 'Omega0', 'alphatau', 'offset')
 