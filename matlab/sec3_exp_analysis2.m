% Computing determinant mean and variances at decreasing lambdas.
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

% Fixed parameters (use the same ones as sec3_exp_analysis)
N = 60;
g = 1.5;
omega0 = 1.0;

par = struct('N', N, 'g', g, 'omega0', omega0);

% Number of instances of det samplings
num_trials = 2000;

% Variable parameters
taum_step = 0.5;
taum_end = 8;
taum_arr = taum_step:taum_step:taum_end;

logz_step = 0.5;
logz_end = 4;
logz_arr = logz_step:logz_step:logz_end;

n_taum = numel(taum_arr);
n_logz = numel(logz_arr);
n_loops = n_taum * n_logz;

% Theoretical sync. frequency
Omega_arr = zeros(n_taum,1);

for n = 1:n_taum
    taum = taum_arr(n);
    mu = 1 / taum;
    p = [1 -omega0 mu^2+g*mu -omega0*mu^2];
    eqOmegas = roots(p).';

    % Keep the real root (if any)
    for k=1:numel(eqOmegas)
        if isreal(eqOmegas(k))
            Omega = real(eqOmegas(k));
        end
    end
    Omega_arr(n) = Omega;
end

% Arrays
pvalues = NaN(n_taum, n_logz);
detMeans = NaN(n_taum, n_logz);
detVars = NaN(n_taum, n_logz);
avgCosVars = NaN(n_taum, n_logz);
avgVarCosExps = NaN(n_taum, n_logz);

% MAIN LOOP:

% Wait bar
f = waitbar(0,'Starting trials...') ;

for i = 1:n_taum
    for j = 1:n_logz
        
        taum = taum_arr(i);
        Omega = Omega_arr(i);
        logz = logz_arr(j);
        z = 10^(-logz);
        
        % Waitbar
        trial_num = (i-1)*n_logz + j;
        waittext = ['trial ' num2str(trial_num) ' out of ' num2str(n_loops)] ;
        waitprog = trial_num / n_loops ;
        waitbar(waitprog, f, waittext) ;
            
        % Implement estimate function
        EST = detEstimatesExp(z, taum, Omega, num_trials, par);
        
        % Fill arrays
        pvalues(i,j) = EST.pvalue;
        detMeans(i,j) = EST.detMean;
        detVars(i,j) = EST.detVar;
        avgCosVars(i,j) = EST.avgCosVar;
        avgVarCosExps(i,j) = EST.avgVarCosExp;
        
    end
end

close(f);

% EXPORT

% Reconfigure arrays
taum = taum_arr.';
logz = logz_arr.';
Omega = Omega_arr.';

% Save
filename = ['detEstimates.mat'] ;
dir_file = fullfile(dir_folder, filename) ;
save(dir_file, 'N', 'g', 'omega0', 'taum', 'logz', 'Omega', 'pvalues', ...
    'detMeans', 'detVars', 'avgCosVars', 'avgVarCosExps');