% Export directory
foldername = 'incoherent_state_dirac' ;
trial = 1; % Increase this

% Generate sample plot
is_fig = true;

% Roots to the incoherent state with dirac distribution

% Parameters
w0 = pi/2;
u_arr = linspace(0, pi, 100);
K = 6;
k_arr = -K:K;

tau_arr = zeros(size(k_arr));
g_arr = zeros(size(k_arr));

ind = 1;
for i=1:numel(u_arr)
    for j=1:numel(k_arr)
        tau = w0^(-1)*(u_arr(i) + pi/2 + 2*pi*k_arr(j));
        g = 2*u_arr(i) / tau;
        
        if (tau > 0) && (g > 0)
            tau_arr(ind) = tau;
            g_arr(ind) = g;
            ind = ind + 1;
        end
        
        if numel(tau_arr) <= ind
            tau_arr = cat(2, tau_arr, zeros(size(k_arr)));
            g_arr = cat(2, g_arr, zeros(size(k_arr)));
        end
    end 
end

% Sort
[~, ord_inds] = sort(tau_arr);
tau_arr = tau_arr(ord_inds);
g_arr = g_arr(ord_inds);

% PLOT

if is_fig
    figure;
    axis on;
    plot(tau_arr, g_arr);
    ylim([0,3]);
    yline(0);
    xline(0);
end