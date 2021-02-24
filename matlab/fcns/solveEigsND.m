function output = solveEigsND(Omega, delta, options)
% Returns a structure with the following elements:
%  - re: real mesh
%  - im: imag mesh
%  - fval: integral mesh
%  - found: array of complex points < tol

% Important options:
%  - rerange: range of real part
%  - imrange: range of imag part
%  - renum: number of steps for real part
%  - imnum: number of steps for imag part
%  - roundingthreshold, numberContours

% Uses the rescaled version of the integral equation z -> gain*z

% Equation parameters
gain = options.gain;
tau0 = options.tau0;
g = options.g;
R = options.scale;

% Numerical parameters
tol = options.tol;
roundingthreshold = options.roundingthreshold;

unum = options.renum;
vnum = options.imnum;
umin = options.rerange(1);
umax = options.rerange(2);
vmin = options.imrange(1);
vmax = options.imrange(2);

% Equation functions
tauE = @(x) tau0/(R*gain) + x/R;

if delta == 0
    F = @(z) g*gain*R*cos(Omega*tau0)*(exp(-z*tau0/(R*gain))-1);
else
    int1 = @(x,z) cos(-Omega*gain*tauE(x) + x).*exp(-z*tauE(x) - x.^2/2/delta^2) / sqrt(2*pi*delta^2);
    int2 = @(x,z) -cos(-Omega*gain*tauE(x) + x).*exp(-x.^2/(2*delta^2))/sqrt(2*pi*delta^2);

    % the non-linear equation is expressed as lambda = F(lambda)
    F = @(z) g*gain*R*integral(@(x) int1(x,z), 0, Inf) + g*gain*R*integral(@(x) int2(x,z), 0, Inf);
end

% evaluate on a grid to find good initial conditions
[u,v] = meshgrid(linspace(umin,umax,unum),linspace(vmin,vmax,vnum));
Fvals = zeros(size(u));

% Wait bar
f = waitbar(0,'Starting trials...') ;
waitk = 0;

for indre = 1:vnum
    for indim = 1:unum
        
        % Waitbar
        waitk = waitk + 1;
        waittext = ['re-iter = ' num2str(indre), ', im-iter = ' num2str(indim)] ;
        waitprog = waitk / (unum * vnum);
        waitbar(waitprog, f, waittext) ;
    
        Fvals(indre, indim) = F(u(indre, indim) + 1i*v(indre, indim));
    end
end

close(f);

% Use root finding for initial conditions meeting the tolerance on the grid
ind = find(log10(abs( u(:) + 1i*v(:) - Fvals(:))) < tol);

% Remove all inds such that u >= 0 and sqrt(u^2 + v^2) < 2*R*g
ind2 = find((u(:) >= 0) & (abs(u(:) + 1i*v(:)) > 2*R*gain*g));
ind = setdiff(ind, ind2);

% Display search point number:
disp(['Number of search points: ' num2str(numel(ind))]);

% disp(numel(ind));
lambda_found = NaN(numel(ind));
residual = NaN(numel(ind));

% Parallel loop
parfor ic=1:numel(ind)
    z_re0 = u(ind(ic));
    z_im0 = v(ind(ic));
    [z,residual(ic)] = fminsearch(@(z) log10(abs(z(1)+1i*z(2) - F(z(1)+1i*z(2)))), [z_re0,z_im0], options);
    lambda_found(ic) = z(1)+1i*z(2);
end

% select out unique solutions
lambda_found = unique(round(lambda_found, roundingthreshold));
lambda_found = lambda_found(~isnan(lambda_found));

residual = unique(round(residual, roundingthreshold));
residual = residual(~isnan(residual));

errors = zeros(size(lambda_found));
for k = 1:numel(lambda_found)
    errors(k) = abs(lambda_found(k) - F(lambda_found(k)));
end
% Store values in output
output = struct();


output.re = u;
output.im = v;
output.fvals = Fvals;
output.found = lambda_found;
output.residual = residual;
output.errors = errors;

end

