function sol = solveDelayKuramoto(parameters, options)
    % Solves the Kuramoto model with a distribution of constant delays.
    % Requires an N by N matrix of constant delays.
    
    % Parameters
    tau = parameters.tau;
    
    N = round(sqrt(numel(tau)));
    w0 = parameters.w0;
    omega = w0*ones(N,1);
    g = parameters.g;
    
    % Starting and end times
    t0 = parameters.t0;
    tf = parameters.tf;

    % Connection topology
    A = g*parameters.A / N;
    
    % Initial condition
    histX = parameters.hist;
    hist_lin = @(t) histX(t-t0) ;
    
    % Functions
    kuraf = @(t,X,Z) modelrhs(t,X,Z,omega,A,N) ;
    tau = tau(:); % Unravel
    
    % solve
    sol = ddesd(kuraf, tau, hist_lin, [t0,tf], options) ;
   
end

function dXdt = modelrhs(t,X,Z,omega,A,N)
    thetadelay = reshape(Z(kron(eye(N),ones(1,N))==1),N,N);
    dXdt = omega + sum( A.*sin( thetadelay - repmat(X,1,N)), 2);
end