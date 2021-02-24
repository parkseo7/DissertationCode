% Set up directory
foldername = 'sec4_2D_numerics4' ; % Change folder name for another trial
cwd = pwd ;
dir_folder = fullfile(cwd, 'data', foldername) ;

if ~exist(dir_folder, 'dir')
   mkdir(dir_folder)
end

% Add to function path
addpath('fcns');
addpath('data');

% Constant parameters
par = struct;
par.w0 = 1.0;
par.g = 1.5/2;
par.t0 = 0;
par.gain = 30;
par.tau0 = 2.0;
par.epsilon = 0.01;
par.offset = 0.0;

% Fixed initial point (change these accordingly)
Omega0 = 0.95;
Delta0 = 0.2;

basetf = 80;

% Variable parameters (dummy values)
par.tf = 50 ;
par.alphatau = 1.0;

% Export names
g = par.g ;
gain = par.gain ;
omega0 = par.w0 ;
tau0 = par.tau0.' ;
alphatau = par.alphatau;
t0 = par.t0 ;
tf = par.tf ;
offset = par.offset;

% Varying alphatau using log scale
alphalog = linspace(0, 2, 10);
alphaarr = 10.^(-alphalog);
tfarr = round(basetf*alphaarr.^(-1), 0);

% DDE options
ddeopts = ddeset() ;
% ddeopts.NormControl = 'on';
ddeopts.MaxStep = 1.0 ;

% Wait bar
f = waitbar(0,'Starting trials...') ;

% MAIN LOOP
for k = 1:numel(alphalog)
    par.alphatau = alphaarr(k);
    par.tf = tfarr(k);
    
    % Waitbar
    waittext = ['alphatau = ' num2str(par.alphatau)] ;
    waitprog = k / numel(alphalog);
    waitbar(waitprog, f, waittext);
    
    % Solve model
    par.hist = IVPhistory2D(Omega0, Delta0, par);
    sol = solvemodel2D(par, ddeopts) ;
    
    % Export (transpose all matrices)
    t = sol.x.' ;
    y = sol.y(1:2,:).' ;
    yp = sol.yp(1:2,:).' ;

    tau = sol.y(3:end,:).' ;
    taup = sol.yp(3:end,:).' ;
    
    alphatau = par.alphatau;
    tf = par.tf;
    
    % Save file
    filename = ['2D_num' num2str(k) '.mat'] ;
    dir_file = fullfile(dir_folder, filename) ;
    save(dir_file, 't', 'y', 'yp', 'tau', 'taup', 'tau0', 'gain', 'omega0', ...
            'g', 'tf', 'Delta0', 'Omega0', 'alphatau', 'offset')
end

close(f)