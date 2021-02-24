% solves the non-linear equation for lambda
clear;
% warning('off','all');

% Parameters
options = struct();
options.gain = 50; % 80
options.tau0 = 0.1;
options.g = 1.5;
options.omega0 = 1.0;

% Delta value
delta = 0.15;
Omega = solveOmega(delta, options);

% Power terms
g = options.g;
tau0 = options.tau0;
gain = options.gain;
tauE = @(x) tau0 + gain*x;
gauss = @(x) exp(-x.^2/(2*delta^2))/sqrt(2*pi*delta^2);
gauss_fun = @(x, m) cos(-Omega*tau0+(1-Omega*gain)*x).*x.^m.*exp(-x.^2/(2*delta^2))/sqrt(2*pi*delta^2);

% Integral of power term
deg_max = 41;
pow_terms = zeros(deg_max+1,1);

for deg = 0:deg_max
    new_int = integral(@(x) gauss_fun(x,deg), 0, inf);
    pow_terms(deg+1) = (-1)^(deg+1)*g*gain*new_int/factorial(deg);
end

% Adjust
pow_terms(1) = 0;
pow_terms(2) = pow_terms(2) + 1;

% Compute and plot poly zeros
poly_roots = roots(flip(pow_terms));

figure;
scatter(real(poly_roots), imag(poly_roots));
yline(0);
xline(0);

