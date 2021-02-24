function phi = IVPhistoryND(freqs, phases, tau0, parameters)
% Returns a history function (with respect to t < 0, with the following
% regions:
% -tau0 < t <= 0: hermite interpolated polynomial
% t < -tau0: linear function with slope freqs
% t = 0: phi(0) = phases
% Here, freqs and phases are 1D row arrays of the same length.

% Parameters
[N,~] = size(phases);
g = parameters.g;
w0 = parameters.w0;
A = parameters.A;

% Define linear functions for t < -tau0
phi0 = @(t) t*freqs + phases;

% Use Kuramoto DDE to obtain derivatives at t = 0:
phi_tau0 = phases - freqs*tau0;
freq = mean(freqs);
dphidt = w0 + (g/N)*sum(A.*sin(bsxfun(@minus, phases, phases')) - freq*tau0,2).';

% Hermite interpolation
tau_min = min(tau0(:));
phi_int = @(t) cubic_int(t, -tau_min, phi_tau0.', 0, phases.', freqs.', dphidt.').';
phi = @(t) hist_pw(t, -tau_min, phi0, phi_int).';

end

function y = hist_pw(t, t0, fun1, fun2)

if t <= t0
    y = fun1(t);
else
    y = fun2(t);
end

end

function [yint,ypint] = cubic_int(tint,t,y,tnew,ynew,yp,ypnew)

h = tnew - t;
s = (tint - t)/h;
s2 = s .* s;
s3 = s .* s2;
slope = (ynew - y)/h;
c = 3*slope - 2*yp - ypnew;
d = yp + ypnew - 2*slope;
yint = y(:,ones(size(tint))) + (h*d*s3 + h*c*s2 + h*yp*s);        
if nargout > 1
  ypint = yp(:,ones(size(tint))) + (3*d*s2 + 2*c*s);  
end   

end
