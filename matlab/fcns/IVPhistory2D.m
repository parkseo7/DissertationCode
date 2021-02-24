function phi = IVPhistory2D(freq, phase, par)
% Returns a history function (with respect to t < 0, with the following
% regions:
% -tau0 < t <= 0: hermite interpolated polynomial
% t < -tau0: linear function with slope freqs
% t = 0: phi(0) = phases
% Here, freq and phase are the desired frequency and phase at t = 0.

tau0 = par.tau0_0;
g = par.g/2;
w0 = par.omega0;
offset = par.offset;

% Use Kuramoto DDE to obtain the phases at t = -tau for the hists.
phi2_tau0 = asin((freq - w0) / g);
phi1_tau0 = asin((freq - w0) / g) + phase;
phi_tau0 = [phi1_tau0, phi2_tau0];
freqs = [freq, freq];

% Define linear functions for t < -tau0
phi0 = @(t) (t - tau0)*freqs + phi_tau0;

% Hermite interpolation
phi_int = @(t) cubic_int(t, -tau0, phi_tau0.', 0, [0, phase].', freqs.', phi_tau0.').';
phi = @(t) hist_pw(t, -tau0, phi0, phi_int).' + offset*[1, 1].';

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

