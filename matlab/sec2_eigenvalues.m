% Clear
clear;

% Add to function path
addpath('data');
addpath('fcns');

% Set up directory (check if it exists)
foldername = 'sec2_eigenvalues' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Parameters
g = -1.5;
tau = 200;
Dplot = [0.12 0.20 0.24 0.40];
Rplot = zeros(size(Dplot));
u0plot = zeros(size(Dplot));

K = 10;
k_arr = -K:K;
eigenvalues = zeros(numel(Dplot), numel(k_arr));

for i = 1:numel(Dplot)
    D = Dplot(i);
    
    % Fixed point:
    errFun = @(u) abs(0.5*g*(1 + erf(u/sqrt(2*D))) - u);
    u0 = fminsearch(errFun, 0);

    % Eigenvalues
    R = g*exp(-u0^2/(2*D))/(sqrt(2*pi*D));
    
    % Solve for eigenvalues
    z_arr = zeros(size(k_arr));
    z_arr2 = zeros(size(k_arr));
    X0 = R*tau*exp(tau);
    for j=1:numel(k_arr)
        indk = k_arr(j);
        z_arr(j) = lambertw(indk, X0)/tau - 1;
    end
    
    % Store
    uplot(i) = u0;
    Rplot(i) = R;
    eigenvalues(i,:) = z_arr;
end


% R-D curve
D_arr = 0.01:0.001:1.0;
R_arr = NaN(size(D_arr));
W0_arr = zeros(size(D_arr));

for i=1:numel(D_arr)
    D = D_arr(i);
    
    % Fixed point:
    errFun = @(u) abs(0.5*g*(1 + erf(u/sqrt(2*D))) - u);
    u0 = fminsearch(errFun, 0);
    
    % Eigenvalues
    R = g*exp(-u0^2/(2*D))/(sqrt(2*pi*D));
    R_arr(i) = R;
    W0_0 = R*tau*exp(tau);
    W0_arr(i) = lambertw(0, W0_0)/tau - 1;
end

% Determine smallest
[~, min_i] = min(abs(W0_arr));
R_min = R_arr(min_i);
D_min = D_arr(min_i);

% Export

% Save
filename = 'eigs2.mat';
dir_file = fullfile(dir_folder, filename);
save(dir_file, 'g', 'tau', 'Dplot', 'Rplot', 'u0plot', 'D_arr', 'R_arr', 'eigenvalues', 'R_min', 'D_min');


