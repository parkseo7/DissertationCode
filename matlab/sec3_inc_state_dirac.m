% Export directory
foldername = 'incoherent_state_dirac' ;
trial = 1; % Increase this

% Generate sample plot
is_fig = true;

% Roots to the incoherent state with dirac distribution

% Parameters
w0 = 1.0;
g = 1.5;
tau0 = 0.5;

% Transcendental parameters
a = -1j*w0;
b = g/2;
c = -tau0;

X0 = -b*c*exp(a*c);

% Solve for eigenvalues
K = 20;
k_arr = -K:K;
z_arr = zeros(size(k_arr));
for i=1:numel(k_arr)
    indk = k_arr(i);
    z_arr(i) = a - (1/c)*lambertw(indk, X0);
end

% EXPORT

% Set up directory (check if it exists)
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Save file
filename = ['sol_' num2str(trial) '.mat'] ;
dir_file = fullfile(dir_folder, filename) ;
save(dir_file, 'z_arr');

% Show preview plots

if is_fig
    figure;
    axis on;
    scatter(real(z_arr), imag(z_arr), '*');
    yline(0);
    xline(0);
end
