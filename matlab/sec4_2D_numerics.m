% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec4_2D_numerics3' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

trialname = 'alpha3';
dir_trial = fullfile(dir_folder, trialname);

if ~exist(dir_trial, 'dir')
   mkdir(dir_trial)
end

% Constant parameters
g = 1.5; % Divided by 2 in the solver.
omega0 = 1.0;
tau0_0 = 0.1;
tau0 = tau0_0*ones(2);
gain = 30;
A = [0 1; 1 0];
alphatau = 0.1; % Toggle this

t0 = 0;
tf = 500;
epsilon = 0.01;
offset = 0.0;

asy = 0.1; % Asymptotic values

% Parameters
parameters = struct('g', g, 'omega0', omega0, 'tau0', tau0, 'tau0_0', tau0_0, ...
    'gain', gain, 'alphatau', alphatau, 't0', t0, 'tf', tf, 'A', A, 'epsilon', epsilon, ...
    'offset', offset, 'phirange', [0,pi], 'numphis', 100);

% Initial conditions grid
step_Omega = 0.1*(g/2);
step_Delta = 0.1;

n_Omega = 21;
L_Delta = 1.0;

Omega_arr = linspace(0,1, n_Omega) - 0.5;
Omega_arr = omega0 + g * Omega_arr;

Delta_arr = step_Delta:step_Delta:L_Delta;
n_Delta = numel(Delta_arr);

total = n_Delta * n_Omega;

% DDE options
options = ddeset() ;
options.OutputFcn = @ddewbar;
% ddeopts.NormControl = 'on';
options.MaxStep = 2.0 ;

% Wait bar
WAIT = waitbar(0,'Starting trials...') ;

% MAIN LOOP
for k = 1:n_Delta
    for l = 1:n_Omega
       
        Delta0 = Delta_arr(k);
        Omega0 = Omega_arr(l);

        % Waitbar
        
        trial_num = (k-1)*n_Omega + l;
        waittext = ['Delta = ' num2str(Delta0) ', Init.freq = ' num2str(Omega0)] ;
        waitprog = trial_num / total;
        waitbar(waitprog, WAIT, waittext) ;

        % Solve model
        freqs0 = Omega0*[1 1];
        phases0 = [0 Delta0];
        histfun = IVPhistory(freqs0, phases0, parameters);
        phifun = histfun.y;
        phipfun = histfun.yp;
        
        % parameters.hist = @(t) freqs0*t + phases0;
        parameters.hist = IVPhistory2D(Omega0, Delta0, parameters); % phifun;
        
        sol = solveDelayKuramoto2D(parameters, options) ;

        % Export (transpose all matrices, format=(index, time))
        t = sol.x.' ;
        y = sol.y(1:2,:).' ;
        yp = sol.yp(1:2,:).' ;
        
        tau = sol.y(3:end,:).' ;
        taup = sol.yp(3:end,:).' ;
        
        % Asymptotic frequency and phases
        t_asy = t(t > (1-asy)*tf);
        y_asy = y(t > (1-asy)*tf,:);
        yp_asy = yp(t > (1-asy)*tf,:);

        t_diff = max(t_asy) - min(t_asy);

        Omega_asy = mean(trapz(t_asy, yp_asy)) / t_diff;
        phases_asy = trapz(t_asy, y_asy - Omega_asy*t_asy) / t_diff;
        Delta_asy = abs(phases_asy(2) - phases_asy(1));
        
        % History arrays
        hist_steps = 1000;
        histt = linspace(t0 - tf/10, t0, hist_steps);
        histy = NaN(2, numel(histt));
        histyp = NaN(2, numel(histt));
        for m=1:numel(histt)
            histy(:,m) = phifun(histt(m));
            histyp(:,m) = phipfun(histt(m));
        end
        histt = histt.';
        histy = histy.';
        histyp = histyp.';
        
        % Save file
        filename = [num2str(k-1, '%02.f') '_' num2str(l-1, '%02.f') '.mat'] ;
        dir_file = fullfile(dir_trial, filename) ;
        save(dir_file, 't', 'y', 'yp', 'tau', 'taup', 'tau0', 'gain', 'omega0', ...
            'g', 't0', 'tf', 'Delta0', 'Omega0', 'alphatau', 'histt', 'histy', ...
            'histyp', 'asy', 'Omega_asy', 'Delta_asy', 'phases_asy');
    end
end

close(WAIT)

% Export parameters
Omega0_arr = Omega_arr.';
Delta0_arr = Delta_arr.';

filename0 = 'parameters.mat';
dir_file0 = fullfile(dir_trial, filename0);
save(dir_file0, 'g', 'omega0', 'tau0', 'gain', 't0', 'tf', 'alphatau', ...
    'epsilon', 'offset', 'Omega0_arr', 'Delta0_arr', 'asy');

