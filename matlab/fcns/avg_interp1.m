function u = avg_interp1(x, v)
% Returns a function structure, with the linear x array and the
% interpolated y-values.

u = struct;
avg_step = 0.5*sum(x(end) - x(1))/size(x,2);

u.x = x(1,1)*ones(1,2*size(x,2)) + avg_step*(1:2*size(x,2)) ;
u.y = interp1(x, v, u.x) ;
u.avg_step = avg_step;

end

