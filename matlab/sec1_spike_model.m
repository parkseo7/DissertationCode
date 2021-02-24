% Section 1: Hudgkin Huxley model
clear;
warning('off','all');

% Add to function path
addpath('fcns');
addpath('data');

% Set up directory (check if it exists)
foldername = 'sec1_HH_model' ;
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% PARAMETERS
Ne = 800; Ni = 200; % Exhitatory, inhibitory neurons

re = rand(Ne, 1); % RS or CH cells
ri = rand(Ni, 1); % uniform sampled from [0,1]

a = [0.02*ones(Ne,1); 0.02+0.08*ri];
b = [0.2*ones(Ne,1); 0.25-0.05*ri];
c = [-65+15*re.^2; -65*ones(Ni,1)];
d = [8-6*re.^2; 2*ones(Ni,1)];

S = [0.5*rand(Ne+Ni,Ne), -rand(Ne+Ni,Ni)];

% Initial values
v0 = -65*ones(Ne+Ni,1);
u0 = b.*v0;

% Time 
T = 100;

% Solution arrays
firings = zeros(Ne+Ni, T);
v = zeros(Ne+Ni, T);
u = zeros(Ne+Ni, T);

% Exclude initial points
v1 = v0;
u1 = u0;

% Wait bar
f = waitbar(0,'Starting trials...') ;

v1_prev = -65;
t_arr = 1:T;
for t = t_arr % simulation of 1000 ms
    
    % Wait progress
    waittext = ['Time: ' num2str(t) 'ms'] ;
    waitprog = t / T;
    waitbar(waitprog, f, waittext) ;
    
    I = [5*randn(Ne,1) ; 2*randn(Ni,1)]; % thalamic input, Gaussian
    fired = find(v1 >= 30);
    firings(:,t) = (v1 >= 30);
    
    v1(fired) = c(fired);
    u1(fired) = u1(fired) + d(fired);
    I = I + sum(S(:,fired), 2);
    
    v1 = 0.5*(0.04*v1.^2 + 5*v1 + 140 - u1 + I);
    v1 = 0.5*(0.04*v1.^2 + 5*v1 + 140 - u1 + I);
    u1 = u1 + a.*(b.*v1 - u1);
    
    % STORE
    v(:,t) = min(v1, 30);
    u(:,t) = u1;
    
end

close(f)


% Save file
filename = 'trial2.mat' ;
dir_file = fullfile(dir_folder, filename);

save(dir_file, 't_arr', 'u', 'v', 'firings')

% figure;
% plot(t_arr, v(1,:));