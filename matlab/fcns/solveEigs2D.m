function output = solveEigs2D(Omega, Delta, parameters, options, fminoptions)
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
gain = parameters.gain;
tau0 = parameters.tau0;
g = parameters.g;

% Numerical parameters
tol = options.tol;
roundingthreshold = options.roundingthreshold;

uR = options.rerange; % real value range
vR = options.imrange; % imag value range

unum = options.renum;
vnum = options.imnum;

% Equation functions
tauE = tau0 + gain*sin(Delta);
C1 = g*cos(-Omega*tauE+Delta);
C2 = g*cos(Delta);
gaintilde = 1 - Omega*gain*cos(Delta);

% Coefficients of polynomials from highest to lowest term
coeffP = [1 (1+C1+C2) (gaintilde*C1+(1+C1)*C2) C1*C2];
coeffQ = [0 0 -C1*C2 -C1*C2];

% Compute polynomial roots
poly_roots = roots(coeffP+coeffQ);

% Eigenvalue function
eig_fun = @(z) sum(coeffP.*[z^3 z^2 z 1] + coeffQ.*[z^3 z^2 z 1]*exp(-z*tauE));

% Evaluate on a grid to find good initial conditions
[u,v] = meshgrid(linspace(-uR,uR,unum+1),linspace(-vR,vR,vnum+1));
Fvals = zeros(size(u));

f = waitbar(0,'Starting trials...') ; % Wait bar
waitk = 0;

for indre = 1:vnum
    for indim = 1:unum
        
        % Waitbar
        waitk = waitk + 1;
        waittext = ['re-iter = ' num2str(indre), ', im-iter = ' num2str(indim)] ;
        waitprog = waitk / (unum * vnum);
        waitbar(waitprog, f, waittext) ;
    
        Fvals(indre, indim) = eig_fun(u(indre, indim) + 1i*v(indre, indim));
    end
end

close(f)

% Use root finding for initial conditions meeting the tolerance on the grid
ind = find(abs(Fvals(:)) < tol);
lambda_found = NaN(numel(ind));
residual = NaN(numel(ind));

% Parallel loop
parfor ic=1:numel(ind)
    z_re0 = u(ind(ic));
    z_im0 = v(ind(ic));
    
    errFun = @(z) log10(abs(eig_fun(z(1)+1i*z(2))));
    
    [z,residual(ic)] = fminsearch(errFun, [z_re0,z_im0], fminoptions);
    lambda_found(ic) = z(1)+1i*z(2);
end

% select out unique solutions
lambda_found = unique(round(lambda_found, roundingthreshold));
lambda_found = lambda_found(~isnan(lambda_found));

residual = unique(round(residual, roundingthreshold));
residual = residual(~isnan(residual));

% Errors
foundVals = NaN(size(lambda_found));
for i=1:numel(lambda_found)
    z0 = lambda_found(i);
    err0 = abs(eig_fun(z0));
    foundVals(i) = err0;
end

% Store values in output
output = struct();

output.re = u;
output.im = v;
output.fvals = Fvals;
output.found = lambda_found;
output.residual = residual;
output.foundVals = foundVals;
output.polyroots = poly_roots;
end



