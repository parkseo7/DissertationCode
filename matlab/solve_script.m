% solves the non-linear equation for lambda
clear;
% warning('off','all');

% Add to function path
addpath('fcns');

% Parameters
options = struct();
options.gain = 50; % 80
options.tau0 = 0.1;
options.g = 1.5;
options.omega0 = 1.0;

% Delta value
delta = 0.015;
Omega = solveOmega(delta, options);

% EIGENVALUES

% Numerical parameters
options.tol = 1.5;
options.roundingthreshold = 1;
options.numberContours = 2^8;
options.MaxIntervalCount = 1e+6;
% options.rerange = 0; % real value range
options.imrange = 200*options.g; % imag value range

options.renum = 60;
options.imnum = 60;
options.remin = -150*options.g;
options.scale = 1;

% Obtain distribution
dist = solveEigsND(Omega, delta, options);
eig_points = dist.found;
eig_res = dist.residual;

actual_eigs = dist.found*options.gain*options.scale;

% plot
figure();
contourf( dist.re, dist.im, log10(1+abs(dist.re + 1i*dist.im - dist.fvals)), options.numberContours, 'EdgeColor', 'none' )
hold on;
plot(real(eig_points),imag(eig_points),'rx')
%plot(lambdar(ind),lambdai(ind),'k.') % initial guesses
xlabel('$$Re \lambda$$','Interpreter','latex')
ylabel('$$Im \lambda$$','Interpreter','latex')
set(gca,'TickLabelInterpreter','Latex');
title('$$log_{10}| \lambda - F(\lambda)|$$','Interpreter','latex')
cbh=colorbar;
set(cbh,'TickLabelInterpreter','Latex')
box on

%figure();
% contourf( dist.re, dist.im, abs(dist.errs), options.numberContours, 'EdgeColor', 'none' )
% hold on;
% plot(real(init_points),imag(init_points),'rx')
% cbh=colorbar;
% box on