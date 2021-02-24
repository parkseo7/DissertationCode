function output = solveOmegaND(delta, options)
% Returns the corresponding Omega solution with given delta value.
% Uses a double integral, while assuming each phase is Gaussian
% distributed.

% Equation parameters
gain = options.gain;
tau0 = options.tau0;
g = options.g;
omega0 = options.omega0;

% Delta values

% Fixed-point solver
delta2 = delta^2;

tauE = @(x,y) max(tau0+gain*sin(x-y), 0);
gauss = @(x) exp(-x.^2/2/delta2) / sqrt(2*pi*delta2);
int_f = @(Omega, x, y) sin(-Omega*tauE(x,y) + (x-y)).*gauss(x).*gauss(y);
F = @(Omega) g*integral2(@(x,y) int_f(Omega, x,y), 0, Inf, 0, inf);
[Omega0, res0] = fminsearch(@(u) abs(u - omega0 - F(u)), omega0);

% Obtain array
u_arr = linspace(omega0-g, omega0+g, 100);
G_arr = zeros(size(u_arr));

if options.get_arr
    G = @(u) u - omega0 - F(u);

    for i = 1:numel(u_arr)
        G_arr(i) = G(u_arr(i));
    end
end

% Output
output.Omega = Omega0;
output.Omega_arr = u_arr;
output.err_arr = G_arr; % L-R error

end

