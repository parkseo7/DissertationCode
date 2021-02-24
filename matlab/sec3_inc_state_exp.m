% Export directory
foldername = 'sec3_incoherent' ;
trial = 1; % Increase this

% Generate sample plot
is_fig = true;

% Roots to the incoherent state with exponential distribution

% Parameters
g = 1.5;
w0 = 1.0;

% Coefficients
taum_arr = linspace(0, 10, 1000);
taum_arr = taum_arr(2:end-1);
mu_arr = 1./taum_arr;

max_real = zeros(size(m_arr));
min_real = zeros(size(m_arr));
eigenvalues = zeros(numel(m_arr), 2);

for k=1:numel(m_arr)
    mu = mu_arr(k);
    p = [1 mu+1j*w0 1j*w0*mu-g*mu/2];
    eigs = roots(p);
    
    % Add largest real part:
    max_real(k) = max(real(eigs));
    min_real(k) = min(real(eigs));
    eigenvalues(k,:) = eigs;
end

% EXPORT
taum = taum_arr.';
max_real = max_real.';
min_real = min_real.';
eigenvalues = eigenvalues.';

% Set up directory (check if it exists)
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Save file
filename = ['incoherent_exp.mat'] ;
dir_file = fullfile(dir_folder, filename) ;
save(dir_file, 'taum', 'min_real', 'max_real', 'eigenvalues');

% Show preview plots

if is_fig
    figure;
    plot(m_arr, log(max_real), '-')
    hold on;
    plot(m_arr, log(abs(min_real)), '.')
end

