% solves the non-linear equation for lambda
clear;

warning('off','all');

% Add to function path
addpath('fcns');

% Parameters
N = 50; %40;
gain = 80;
tau0 = 0.1;
g = 1.5;
omega0 = 1.0;

options = struct;
% Add parameters
options.gain = gain;
options.tau0 = tau0;
options.g = g;
options.omega0 = omega0;
options.N = N;

% Obtain Omega(delta) array
deltaarr = 0:0.005:0.1;

Omegaarr = omega0*ones(size(deltaarr));
for k = 1:numel(deltaarr)
    Omegaarr(k) = solveOmega(deltaarr(k), options);
end

% Define M_ij term
tauE =  @(x) (tau0 + gain*x).*double((tau0 + gain*x) > 0);
C_ij = @(x, u) cos(-u*tauE(x) + x);
M_ij = @(x, u, z) C_ij(x, u).*((z+1)*exp(-z*tauE(x)) - u*gain*cos(x).*double(tauE(x) > 0));

meanarr = zeros(numel(deltaarr), 1);
stdarr = zeros(numel(deltaarr), 1);

for k = 1:numel(deltaarr)
    z0 = 0.0;
    delta = deltaarr(k);
    Omega = Omegaarr(k);

    EX = integral(@(x) M_ij(x,Omega,z0).*exp(-x.^2/2/delta^2)/sqrt(2*pi), -inf, inf);
    EX2 = integral(@(x) M_ij(x,Omega,z0).^2.*exp(-x.^2/2/delta^2)/sqrt(2*pi), -inf, inf);
    stdX = sqrt(EX2 - EX^2);
    
    meanarr(k) = EX;
    stdarr(k) = stdX;
end

figure
hold on;
plot(deltaarr, meanarr)
plot(deltaarr, stdarr)