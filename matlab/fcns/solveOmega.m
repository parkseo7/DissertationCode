function output = solveOmega(delta, options)
% Returns the corresponding Omega solution with given delta value.

% Equation parameters
gain = options.gain;
tau0 = options.tau0;
g = options.g;
omega0 = options.omega0;

% Delta values

% Fixed-point solver
delta2 = delta^2;

if (delta == 0)
    F = @(Omega) sin(-Omega*tau0);
else
    tauE = @(Delta) max(tau0+gain*Delta, 0);
    gauss = @(Delta) exp(-Delta.^2/2/delta2) / sqrt(2*pi*delta2);
    int_f = @(Omega, Delta) sin(-Omega*tauE(Delta) + Delta).*gauss(Delta);
    F = @(Omega) g*integral(@(Delta) int_f(Omega, Delta), -Inf, Inf);
end

[Omega0, res0] = fminsearch(@(u) abs(u - omega0 - F(u)), omega0);
output = Omega0;

end

