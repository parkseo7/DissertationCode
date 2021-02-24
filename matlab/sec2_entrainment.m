% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec2_entrainment' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

is_plot = false;

% Parameters
g = -1.5;
D = 0.22;
tau = 200;

% Acquire Hopf frequencies using LambertW

% Fixed point:
errFun = @(u) abs(0.5*g*(1 + erf(u/sqrt(2*D))) - u);
u0 = fminsearch(errFun, 0);

% Eigenvalues
R = g*exp(-u0^2/(2*D))/(sqrt(2*pi*D));

% Solve for eigenvalues
K = 10;
k_arr = 0:K;
z_arr = zeros(size(k_arr));
z_arr2 = zeros(size(k_arr));
X0 = R*tau*exp(tau);
for i=1:numel(k_arr)
    indk = k_arr(i);
    z_arr(i) = lambertw(indk, X0)/tau - 1;
    z_arr2(i) = lambertw(-indk, X0)/tau - 1;
end

% Stimulation frequencies and amplitudes
omega0 = 2*pi*10/1000;
omega = abs(imag([z_arr(2) z_arr(3)]));
S = [0.5 0.5];

% Times
t0 = 0;
tf = 6000;
t1 = 2000;
t2 = 4000;

% Model
noiseMolFun = noiseMollifier(t1, t2);
Sfun = @(t) noiseMolFun(t)*sum(S.*cos(omega*t));
u_fun = @(t, u, utau) -u + 0.5*g*(1 + erf(utau/sqrt(2*D))) + Sfun(t);
hist_fun = @(t) 1;

% Solve DDE:
sol = dde23(u_fun, tau, hist_fun, [t0, tf]);

% Generate arrays
Sarr_x = linspace(t0, tf, 1000);

Sarr_y = zeros(size(Sarr_x));
for j=1:numel(Sarr_x)
    Sarr_y(j) = Sfun(Sarr_x(j));
end

% FOURIER TRANFORM
M = 2*numel(sol.x);
u_PreArr = PowerSpectrum(sol, [t0+tau, t1], M);
u_EntArr = PowerSpectrum(sol, [t1+tau, t2-tau], M);
u_PostArr = PowerSpectrum(sol, [t2+tau, tf-tau], M);

logPre = log10(1 + abs(u_PreArr.y));
logEnt = log10(1 + abs(u_EntArr.y));
logPost = log10(1 + abs(u_PostArr.y));

if is_plot
    % PLOT
    figure;
    plot(sol.x, sol.y);

    % PLOT
    figure;
    hold on;
    plot(1000*u_PreArr.x(2:end), logPre(2:end));
    plot(1000*u_EntArr.x(2:end), logEnt(2:end));
    plot(1000*u_PostArr.x(2:end), logPost(2:end));
    xline(1000*omega/(2*pi));
    xlim([0,100]);

    figure;
    hold on;
    scatter(1000*real(z_arr), 1000*imag(z_arr));
    scatter(1000*real(z_arr2), 1000*imag(z_arr2));
    xline(0)
end

% EXPORT

% Compile
t = sol.x;
y = sol.y;
yp = sol.yp;

freqPre = u_PreArr.x;
powPre = u_PreArr.y;
freqEnt = u_EntArr.x;
powEnt = u_EntArr.y;
freqPost = u_PostArr.x;
powPost = u_PostArr.y;

eigs1 = z_arr;
eigs2 = z_arr2;

% Save
filename = 'trial1.mat';
dir_file = fullfile(dir_folder, filename);
save(dir_file, 'g', 'D', 'S', 'tau', 'omega', 'omega0', 'R', 'u0', ...
    't0', 't1', 't2', 'tf', 't', 'y', 'yp', 'Sarr_x', 'Sarr_y', ...
    'eigs1', 'eigs2', 'freqPre', 'powPre', 'freqEnt', 'powEnt', ...
    'freqPost', 'powPost')






