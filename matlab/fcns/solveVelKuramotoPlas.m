function sol = solveVelKuramotoPlas(parameters, options)

    % Parameters
    N = parameters.N ;
    w0 = parameters.w0 ;
    omega = w0*ones(N,1) ;
    g = parameters.g ;
    dist0 = parameters.dist; % Matrix of sampled distances
    distm = parameters.distm;
    vel0_0 = parameters.vel0 ;
    vel_min = parameters.vel_min;
    vel_max = parameters.vel_max;
    
    alphar = parameters.alphatau ;
    betar = parameters.betatau;
    t0 = parameters.t0 ;
    tf = parameters.tf ;
    gain = parameters.gain;
    % ep0 = parameters.epsilon;
    
    A = (g/N)*parameters.A;
    vel0 = vel0_0*ones(N,N);
    % tau0 = dist / vel0_0;
    
    dist_norm = dist0 / distm; % max(dist0(:));
    % initial condition (constantly distributed around half-circle at t0)
    histX = parameters.hist;
    % hist_linX = @(t) offset*(pi/N)*(0:N-1) + omega.'*(t - t0) ;
    % hist_linX = @(t) (pi/N)*(0:N-1) + omega.'*(t - t0) ;
    hist_lin = @(t) packX(histX(t-t0), vel0) ;
    
    % Functions
    kuraf = @(t,X,Z) modelrhs(t,X,Z,omega,A,alphar,betar,gain,vel0_0,dist_norm) ;
    tauf = @(t,X) delays(t,X,dist0,vel_min,vel_max);
    
    % solve
    sol = ddesd(kuraf, tauf, hist_lin, [t0,tf], options) ;
   
end

function X = packX( theta, vel )
    X = [ theta(:) ; vel(:) ];
end

function [theta, vel, N] = unpackX( X )
    N = round( 0.5*(sqrt(4*numel(X)+1)-1) );
    theta = X(1:N);
    vel = reshape( X(N+1:end), [N,N] );
end

function dXdt = modelrhs(t,X,Z,omega,A,alphar,betar,gain,vel0,dist_norm)
    [theta, vel, N] = unpackX( X );
    thetadelay = Z(1:N,:); % Cut-off delay component (theta is column vector)
    thetadelay = reshape(thetadelay(kron(eye(N),ones(1,N))==1),N,N);
    dthetadt = omega + sum( A.*sin( thetadelay - repmat(theta,1,N)), 2);
    % ordert = abs(sum(exp(1j*theta)))/N;
    % meantheta = repmat(sum(sin(bsxfun(@minus,theta',theta)),2),1,N)/N;
    Delta = bsxfun(@minus,theta',theta);
    % ind_bound = isbound(vel, vel_min, vel_max);
    dveldt = alphar*(-betar*dist_norm.*(vel - vel0) + gain*max(sin(Delta),0));
        % + betar*meantheta);
    dXdt = packX( dthetadt, dveldt );
    
end

function d = delays(t, X, dist1,vel_min,vel_max)
    [~, vel] = unpackX( X );
    
    vel = min(max(vel,vel_min),vel_max); % Bound velocity
    tau = dist1 ./ vel;
    d = t - tau(:);
    
    % tau = dist ./ vel;
    disp(['Time = ' num2str(t) ', min = ' num2str(min(tau(:))), ', max = ' num2str(max(tau(:)))]);
    
end

function u = isbound(vel, vel_min, vel_max)

% Adjust rates for tau extremely close to 0, so that they move slower
u = double(vel > vel_min).*double(vel < vel_max);

end
