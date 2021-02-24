function output = findlambda(init, Omega, delta, options)
% Uses fmin to find the exact locations of eigenvalues. init is the
% starting guesses. Returns the following:
% - found: An array of fmin solutions
% - residual: An array of residual (actual fmin values, should be close to
% 0).

% Equation parameters
R = options.scale;
gain = options.gain / R;
tau0 = options.tau0 / R;
g = R * options.g;
delta2 = 2 * delta^2;

% Numerical parameters
tol = options.tol;
roundingthreshold = options.roundingthreshold;
MIC = options.MaxIntervalCount;

uR = options.rerange; % real value range
umin = options.remin; % minimum real value
vR = options.imrange; % imag value range

unum = options.renum;
vnum = options.imnum;

% Equation functions
tauE = @(Delta) max(tau0+gain*Delta, 0);

% fmin arrays
found = NaN(size(init));
residual = NaN(size(init));

% Find values using parallel loops
parfor k=1:numel(init)
    z_re0 = real(init(k));
    z_im0 = imag(init(k));
    INT = @(x,z) g * cos(-Omega*tauE(x) + x).*(exp(-z*tauE(x)-x.^2/delta2) - exp(-x.^2/delta2))/sqrt(pi*delta2);
    % F = @(z) quadgk(@(x) INT(x,z), 0, Inf, 'MaxIntervalCount', MIC);
    F = @(z) integral(@(x) INT(x,z), 0, Inf);
    [z,residual(k)] = fminsearch(@(z) abs(z(1)+1i*z(2) - F(z(1)+1i*z(2))), [z_re0,z_im0] );
    found(k) = z(1)+1i*z(2);
end

% select out unique solutions
found = unique(round(found, roundingthreshold));
found = found(~isnan(found));

residual = unique(round(residual, roundingthreshold));
residual = residual(~isnan(residual));

output = struct();
output.found = found;
output.residual = residual;

end

