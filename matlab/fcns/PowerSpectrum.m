function u = PowerSpectrum(sol, tint, M)
% Returns the power spectrum of a periodic solution from time t = t1 to t =
% t2 using Fast Fourier Series. If time is in ms, multiply x (frequency) by
% 1000.

ts = tint(1);
tf = tint(2);

[~, ts_ind] = min(abs(sol.x - ts));
[~, tf_ind] = min(abs(sol.x - tf));
t_arr = sol.x(ts_ind:tf_ind) - ts;
y_arr = sol.y(ts_ind:tf_ind);


% Interpolate
U = avg_interp1(t_arr, y_arr, M);
U_fft = abs(fft(U.y));
U_N = numel(U.x);
U_L = U.x(end);
post_fft_x = (1:U_N)/U_L;

u.x = post_fft_x;
u.y = U_fft;
u.N = U_N;
u.L = U_L;

end

function u = avg_interp1(x, v, M)
% Returns a function structure, with the linear x array and the
% interpolated y-values.

u = struct;
% avg_step = 0.5*sum(x(end) - x(1))/numel(x);
u.x = linspace(x(1), x(end), M);
% u.x = x(1)*ones(1,2*numel(x)) + avg_step*(1:2*numel(x)) ;
u.y = interp1(x, v, u.x) ;
% u.avg_step = avg_step;

end

