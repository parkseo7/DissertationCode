function output = solveEigsDet(Omega, options)
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

% This is the non-delay version.

% Equation parameters
tau0 = options.tau0; % tau matrix
g = options.g;
taum = options.taum;
mu = 1 / taum;

% Numerical parameters
tol = options.tol;
roundingthreshold = options.roundingthreshold;

umin = options.rerange(1);
umax = options.rerange(2);
vmin = options.imrange(1);
vmax = options.imrange(2);

unum = options.renum;
vnum = options.imnum;

% Define matrices
N = round(sqrt(numel(tau0)));

% Determinant function
detlambda = @(z) det(eye(N) - g*cos(Omega*tau0).*exp(-z*tau0)/N./(z + g*repmat(sum(cos(Omega*tau0),2),1,N)/N));
% detlambda = @(z) det(eye(N) - g*cos(Omega*tau0).*exp(-z*tau0)/N)./(z + g*mu^2/(Omega^2+mu^2));

% evaluate on a grid to find good initial conditions
[u,v] = meshgrid(linspace(umin, umax,unum+1),linspace(vmin, vmax,vnum+1));
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
    
        Fvals(indre, indim) = abs(detlambda(u(indre, indim) + 1i*v(indre, indim)));
    end
end

close(f);

% Use root finding for initial conditions meeting the tolerance on the grid
ind = find(log10(Fvals(:)) < tol);

% Display search point number:
disp(['Number of search points: ' num2str(numel(ind))]);

% disp(numel(ind));
lambda_found = NaN(numel(ind));
residual = NaN(numel(ind));

% Parallel loop
parfor ic=1:numel(ind)
    z_re0 = u(ind(ic));
    z_im0 = v(ind(ic));
    [z,residual(ic)] = fminsearch(@(z) log10(abs(detlambda(z(1)+1i*z(2)))), [z_re0,z_im0], options);
    lambda_found(ic) = z(1)+1i*z(2);
end

% select out unique solutions
lambda_found = unique(round(lambda_found, roundingthreshold));
lambda_found = lambda_found(~isnan(lambda_found));

residual = unique(round(residual, roundingthreshold));

errors = zeros(size(lambda_found));
for k = 1:numel(lambda_found)
    errors(k) = abs(detlambda(lambda_found(k)));
end

lambda_found = lambda_found(~isinf(errors));
errors = errors(~isinf(errors));

cap = options.cap;
lambda_found = lambda_found(errors < cap);
errors = errors(errors < cap);

% Store values in output
output = struct();

output.re = u;
output.im = v;
output.fvals = Fvals;
output.found = lambda_found;
output.residual = residual;
output.errors = errors;
output.zero = abs(detlambda(0));

end