function sol = solveDelayKuramotoInj(parameters, options)

    % Parameters
    N = parameters.N ;
    w0 = parameters.w0 ;
    omega = w0*ones(N,1) ;
    g = parameters.g ;
    tau0_0 = parameters.tau0 ;
    kappa = parameters.gain ;
    alphar = parameters.alphatau ;
    epsilon = parameters.epsilon ;
    t0 = parameters.t0 ;
    tf = parameters.tf ;
    
    tau0 = tau0_0*ones(N,N);
    
    inj = parameters.inj;
    tinjs = parameters.tinjs;
    tinjf = parameters.tinjf;
   
    % Connection matrix
    A_s = (g/N)*ones(N);
    A_f = (g/N)*double(rand(N,N) < inj);
    
    % initial condition (constantly distributed around half-circle at t0)
    histX = parameters.hist;
    % hist_linX = @(t) offset*(pi/N)*(0:N-1) + omega.'*(t - t0) ;
    % hist_linX = @(t) (pi/N)*(0:N-1) + omega.'*(t - t0) ;
    hist_lin = @(t) packX(histX(t-t0), tau0) ;
    
    % Functions
    kuraf = @(t,X,Z) modelrhs(t,X,Z,omega,A_s,A_f,tinjs,tinjf,kappa,alphar,tau0,epsilon) ;
    tauf = @delays;
    
    % solve
    sol = ddesd(kuraf, tauf, hist_lin, [t0,tf], options) ;
    
    % Solution components
    sol.tau0 = tau0 ;
    sol.A_inj = (N/g)*(A_s - A_f); % Which ones are not injured
    
    sol.molarr = zeros(numel(sol.x),1);
    for k = 1:numel(sol.x)
        X = sol.x(k);
        sol.molarr(k) = 1 - mollifier(X, tinjs, tinjf);
    end
end

function X = packX( theta, tau )
    X = [ theta(:) ; tau(:) ];
end

function [theta, tau, N] = unpackX( X )
    N = round( 0.5*(sqrt(4*numel(X)+1)-1) );
    theta = X(1:N);
    tau = reshape( X(N+1:end), [N,N] );
end

function dXdt = modelrhs(t,X,Z,omega,A1,A2,t1,t2,kappa,alphar,tau0,epsilon)
    [theta, tau, N] = unpackX( X );
    thetadelay = Z(1:N,:); % Cut-off delay component
    thetadelay = reshape(thetadelay(kron(eye(N),ones(1,N))==1),N,N);
    A_t = A1 - mollifier(t,t1,t2)*A2;
    dthetadt = omega + sum( A_t.*sin( thetadelay - repmat(theta,1,N)), 2);
    dtaudt = alphar*posind(tau, epsilon).*( -(tau - tau0) + kappa*sin(bsxfun(@minus,theta',theta)));
    dXdt = packX( dthetadt, dtaudt );
    
    disp(['Time = ' num2str(t) ', min = ' num2str(min(tau, [], 'all')), ', max = ' num2str(max(tau, [], 'all'))])
end

function d = delays(t, X)
    [~, tau] = unpackX( X );
    tau_pos = (tau > 0) .* tau; % Ensure positivity
    d = t - tau_pos(:);
end

function u = posind(tau, epsilon)
% Returns a vector with each component being 1 if tau_j > 0 and 0
% otherwise

u = double(tau >= epsilon);

% For each index in x_inds, replace x_pos:
MOL = @(t) exp(-(t-1).^-2).*exp(-(t+1).^-2);
MOL2 = @(v) MOL(-1 + 2*v/epsilon);
INT = integral(@(s) MOL2(s), 0, epsilon);

% With 0 < tau < epsilon, choose the average value.
tau_inds1 = find((tau > 0).*(tau < epsilon));

if numel(tau_inds1) > 0
    tau_avg = mean(tau(tau_inds1));
    u_avg = integral(@(s) MOL2(s), 0, tau_avg) / INT;
    
    u(tau_inds1) = u_avg ;
end
end

function u = mollifier(t, t1, t2)
% Returns a smooth decay value starting at t_injs and ending at t_injf.

MOL = @(s) exp(-(s-1).^-2).*exp(-(s+1).^-2);
MOL2 = @(v) MOL(-1 + 2*v/(t2 - t1));
INT = integral(@(s) MOL2(s), 0, (t2 - t1));

if t < t1
    u = 0;
elseif t > t2
    u = 1;
else
    u = integral(@(s) MOL2(s), 0, t-t1) / INT;
end

end