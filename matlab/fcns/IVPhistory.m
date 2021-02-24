function histfun = IVPhistory(freqs, phases, parameters)
% Returns a history function (with respect to t < 0, with the following
% regions:
% -tau0 < t <= 0: hermite interpolated polynomial
% t < -tau0: linear function with slope freqs
% t = 0: phi(0) = phases
% Here, freqs and phases are 1D row arrays of the same length.

% Parameters
N = numel(phases);
tau0 = parameters.tau0;
g = parameters.g;
w0 = parameters.w0;
A = parameters.A;
t0 = parameters.t0;

% First starting point (min_j tau_ij for each i)
t1 = -min(tau0(:));
p1_arr = phases;
m1_arr = freqs;

% Define velocity at t = 0
t2 = 0;
vel = @(phi, i) w0 + (g/N)*sum(A(i,:).*sin(freqs(i)*(-tau0(i,:)-t1) - phi));
tint = @(t) (t - t1) / (t2 - t1);

% Options for optimization search
phirange = parameters.phirange;
numphis = parameters.numphis;
phi_arr = linspace(phirange(1), phirange(2), numphis);

% Optimize over each index
herm_coeffs = zeros(N, 4);

for i=1:N
    
    % Position and velocity at t = t1
    p1 = p1_arr(i);
    m1 = m1_arr(i);
        
    % Linear coefficients (without interpolation)
    b1 = freqs(i)*(t2 - t1);
    b0 = p1;
    
    linear = @(x) b1*x + b0;
    linear_int = @(t) linear(tint(t));
    
    dist_arr = NaN(size(phi_arr));
    dist_arr2 = NaN(size(phi_arr));
    for k = 1:numel(phi_arr)
        
        % Position and velocity at t = 0
        p2 = phi_arr(k);
        m2 = vel(p2, i);
        
        % Cubic coefficients
        coeffs = hermite(t1, t2, p1, p2, m1, m2);
        a3 = coeffs(1);
        a2 = coeffs(2);
        a1 = coeffs(3);
        a0 = coeffs(4);
        
        % Minimize quadratic:
        a = 3*a3;
        b = 2*a2;
        c = a1 - b1;
        
        % Min point:
        t_c = a2 / 3*a3;
        if (t_c - t1) / (t2 - t1) > 1 || (t_c - t1) / (t2 - t1) < 0
            dist_c = 0;
        else
            dist_c = abs(a*t_c^2 + b*t_c + c);
        end
            
        dist_arr2(k) = max(dist_c, abs(m2 - m1));
        
        cubic = @(x) a3*x^3 + a2*x^2 + a1*x + a0;
        cubic_int = @(t) cubic(tint(t));
        
        % Discriminent (of derivative)
        a = 3*a3;
        b = 2*a2;
        c = a1 - b1;
        disc = b^2 - 4*a*c;
        
        % Find maximum distance
        dist0 = abs(cubic_int(t2) - linear_int(t2));
        
        % No critical points
        if disc <= 0
            dist_arr(k) = dist0;
        
        % Two critical points
        else
            t_c1 = (-b + sqrt(disc)) /2/a;
            t_c2 = (-b - sqrt(disc)) /2/a;
            is_within1 = (-t1 < t_c1 && t_c1 < t2);
            is_within2 = (-t1 < t_c2 && t_c2 < t2);
            
            if ~is_within1 && ~is_within2
                dist_arr(k) = dist0;
            else
                dist1 = abs(cubic_int(t_c1) - linear_int(t_c1));
                dist2 = abs(cubic_int(t_c2) - linear_int(t_c12));
                dist_arr(k) = max(dist0, dist1*is_within1, dist2*is_within2);
            end
            
        end
        
    end
    
    % Use smallest distance index
    [~,ind] = min(dist_arr2);
    
    % Position and velocity at t = 0
    p2 = phi_arr(ind);
    m2 = vel(p2, i);
    
    % Cubic coefficients
    herm_coeffs(i,:) = hermite(t1, t2, p1, p2, m1, m2);
    
end

% Derivative interpolation
hermp_coeffs = zeros(N, 3);
hermp_coeffs(:,1) = 3*herm_coeffs(:,1);
hermp_coeffs(:,2) = 2*herm_coeffs(:,2);
hermp_coeffs(:,3) = herm_coeffs(:,3);

% Define linear functions for t < -t1
linear_fun = @(t) freqs*(t - t1) + phases;
const_fun = @(t) freqs;

% Define cubic interpolation function
cubic_fun0 = @(t) sum(herm_coeffs.*repmat([t^3 t^2 t 1], N,1),2);
cubic_fun = @(t) cubic_fun0(tint(t));

cubicp_fun0 = @(t) sum(hermp_coeffs.*repmat([t^2 t 1], N,1),2);
cubicp_fun = @(t) cubicp_fun0(tint(t)) / (t2 - t1);

phifun = @(t) hist_pw(t-t0, t1-t0, linear_fun, cubic_fun).';
phipfun = @(t) hist_pw(t-t0, t1-t0, const_fun, cubicp_fun).';

histfun.y = phifun;
histfun.yp = phipfun;

end


function y = hist_pw(t, t_int, fun1, fun2)

if t <= t_int
    y = fun1(t);
else
    y = fun2(t);
end

end


function coeffs = hermite(t1, t2, p1, p2, m1, m2)

tm1 = (t2 - t1)*m1;
tm2 = (t2 - t1)*m2;

a3 = 2*p1 + tm1 - 2*p2 + tm2;
a2 = -3*p1 - 2*tm1 + 3*p2 - tm2;
a1 = tm1;
a0 = p1;
    
coeffs = [a3 a2 a1 a0];

end

