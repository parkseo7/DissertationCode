function output = solveOmega2D(parameters, options)
% Find the synchronization frequency (Omega) and corresponding positive
% phase difference of the 2D Kuramoto adaptive model
% The output has the following components:
% - Omega: Synchronization frequencies
% - Delta: Corresponding phase differences
% - err: Corresponding error of the frequencies
% - x_arr: x-array (Omega) for the error function
% - y_arr: y-array (err) for the error function

% Parameters
g = parameters.g;
w0 = parameters.omega0;
tau0 = parameters.tau0;
gain = parameters.gain;

% Options
M = options.numsteps;
roundingthreshold = options.roundingthreshold;

% Error function
err_fun = @(u) (u-w0) - g*sin(-u.*(tau0 + gain*(w0-u)/g) + asin((w0-u)/g));
err_fun_cc = @(u) (u-w0) - g*sin(-u.*(tau0 + gain*(w0-u)/g) + asin(ccat((w0-u)/g)));

% Error function values
x_arr = linspace(w0-g, w0+g, M);
y_arr = err_fun(x_arr);

% Use root finding for initial conditions meeting the tolerance on the grid
ind = find(y_arr(1:end-1).*y_arr(2:end)<0);
found = NaN(numel(ind));

% Use fmin search on [w0-g, w0+g]
parfor k=1:numel(ind)
    [found(k),~] = fzero(err_fun_cc, x_arr(ind(k)));
end

% select out unique solutions
found = unique(round(found, roundingthreshold));
found = found(~isnan(found));

% Store
output = struct();
output.Omega = found;
output.Delta = asin((w0-found)/g);
output.error = abs(err_fun(found));
output.x_arr = x_arr;
output.y_arr = y_arr;

end

function y = ccat(u)
    if u < -1
        y = -1;
    elseif u > 1
        y = 1;
    else
        y = u;
    end
end
