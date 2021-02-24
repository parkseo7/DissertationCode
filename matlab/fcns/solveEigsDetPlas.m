function output = solveEigsDetPlas(Omega, Delta, options)
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

% Equation parameters
N = options.N;
g = options.g;
w0 = options.omega0;
gain =  options.gain;
tau0 = options.tau0;

% Numerical parameters
R = options.scale;
logtol = options.logtol;
roundingthreshold = options.roundingthreshold;

umin = options.rerange(1);
umax = options.rerange(2);
vmin = options.imrange(1);
vmax = options.imrange(2);

unum = options.renum;
vnum = options.imnum;

% Dummy Gaussian values
% phi = abs(normrnd(0, delta*sqrt(1-2/pi), N, 1)); % Normal -> Folded normal std
% Delta = bsxfun(@minus, phi, phi');

% Delta = normrnd(0, delta, N);

tauE = tau0 + gain*sin(Delta);
tauE = tauE.*double(tauE > 0);
HSDelta = double(tauE > 0);

costauE = cos(-Omega*tauE + Delta);

% Determinent function
detlambda = @(z) det(eye(N) - (g/N)*costauE.*((z+1)*exp(-z*tauE) - ...
    Omega*gain*cos(Delta).*HSDelta) ./ ...
    (z*(z+1) + (g/N)*repmat(sum(costauE.*(z + 1 - ...
    Omega*gain*cos(Delta).*HSDelta),2),1,N)));

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
ind = find(log10(Fvals(:)) < logtol);

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

% Select out unique solutions
lambda_found = unique(round(lambda_found, roundingthreshold));
lambda_found = lambda_found(~isnan(lambda_found));

residual = unique(round(residual, roundingthreshold));

% Errors at located points
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
output.Delta = Delta;

end

